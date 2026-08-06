"""Effect-size and context-dependence statistics.

Two independent readouts, deliberately:
  1. Model-based  -- pseudobulk DE log2FC signatures (pseudobulk.py + pydeseq2)
  2. Model-free   -- E-distance with a permutation null (pertpy)

If a perturbation is called a hit by one and not the other, that disagreement
is a finding about the method, and worth a paragraph.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "signature_matrix",
    "cross_condition_correlation",
    "call_context_dependent",
    "guide_enrichment",
]


def signature_matrix(de_results: dict[tuple[str, str], pd.DataFrame],
                     value_col: str = "log2FoldChange",
                     genes: list[str] | None = None) -> pd.DataFrame:
    """Assemble per-contrast DE tables into a (perturbation x condition) x gene matrix.

    de_results keys are (perturbation, condition) tuples. Index of the returned
    frame is a MultiIndex on those. This is the same object you built for the
    Norman genetic-interaction manifold -- only the second axis has changed
    from 'other gene' to 'environment'.
    """
    frames = {}
    for (pert, cond), df in de_results.items():
        frames[(pert, cond)] = df[value_col]
    mat = pd.DataFrame(frames).T
    mat.index = pd.MultiIndex.from_tuples(mat.index, names=["perturbation", "condition"])
    if genes is not None:
        mat = mat.reindex(columns=genes)
    return mat.fillna(0.0)


def cross_condition_correlation(sig: pd.DataFrame, reference: str) -> pd.DataFrame:
    """For each perturbation, correlate its signature in each condition vs reference.

    Low correlation = the perturbation does something DIFFERENT depending on
    environment. That is the primary screen readout of this project.
    """
    out = []
    for pert in sig.index.get_level_values("perturbation").unique():
        block = sig.loc[pert]
        if reference not in block.index:
            continue
        ref = block.loc[reference]
        for cond in block.index:
            if cond == reference:
                continue
            r = np.corrcoef(ref.values, block.loc[cond].values)[0, 1]
            out.append(
                {
                    "perturbation": pert,
                    "condition": cond,
                    "pearson_vs_reference": r,
                    "l2_reference": float(np.linalg.norm(ref.values)),
                    "l2_condition": float(np.linalg.norm(block.loc[cond].values)),
                }
            )
    df = pd.DataFrame(out)
    # >1 means the perturbation's effect is AMPLIFIED in this environment.
    df["magnitude_ratio"] = df["l2_condition"] / df["l2_reference"].replace(0, np.nan)
    return df


def call_context_dependent(corr_df: pd.DataFrame,
                           edist_df: pd.DataFrame,
                           cfg: dict[str, Any]) -> pd.DataFrame:
    """Join the correlation and E-distance readouts and apply config thresholds.

    A perturbation is context-dependent only if it (a) has a REAL effect in at
    least one condition, by permutation-tested E-distance, and (b) that effect
    diverges between conditions. Requiring (a) is what stops noise from
    dominating the hit list -- two uncorrelated null signatures also have low
    correlation.
    """
    t = cfg["context"]["call_thresholds"]
    merged = corr_df.merge(edist_df, on=["perturbation", "condition"], how="left")
    merged["real_effect"] = merged["edist_pvalue"] < t["min_edistance_significance"]
    merged["diverges"] = merged["pearson_vs_reference"] < t["max_cross_condition_pearson"]
    merged["context_dependent"] = merged["real_effect"] & merged["diverges"]
    return merged.sort_values("pearson_vs_reference")


def guide_enrichment(adata, cfg: dict[str, Any]) -> pd.DataFrame:
    """Perturbation abundance shift between each condition and the reference.

    In the TIL co-culture arm cells are under selection, so a perturbation's
    representation is itself a phenotype: immune-evasion hits should be
    enriched, sensitisers depleted. This is a classic pooled-screen readout
    that costs nothing extra given the data is already loaded -- and it
    answers a different question from the expression readout.
    """
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests

    s = cfg["schema"]["obs"]
    ref = cfg["schema"]["conditions"]["reference"]
    pc = int(cfg["selection"]["pseudocount"])

    tab = (
        adata.obs.groupby([s["perturbation"], s["condition"]], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    totals = tab.sum(axis=0)

    rows = []
    for cond in tab.columns:
        if cond == ref:
            continue
        for pert in tab.index:
            a, b = tab.loc[pert, cond] + pc, totals[cond] - tab.loc[pert, cond] + pc
            c, d = tab.loc[pert, ref] + pc, totals[ref] - tab.loc[pert, ref] + pc
            odds, p = fisher_exact([[a, b], [c, d]])
            rows.append(
                {
                    "perturbation": pert,
                    "condition": cond,
                    "n_condition": int(tab.loc[pert, cond]),
                    "n_reference": int(tab.loc[pert, ref]),
                    "log2_odds_ratio": float(np.log2(odds)) if odds > 0 else np.nan,
                    "pvalue": p,
                }
            )
    df = pd.DataFrame(rows)
    df["padj"] = multipletests(
        df["pvalue"], method=cfg["selection"]["fdr_method"]
    )[1]
    return df.sort_values("log2_odds_ratio")
