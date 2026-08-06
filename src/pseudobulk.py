"""Pseudobulk aggregation for condition-stratified differential expression.

The design constraint that governs this whole module: every contrast is
perturbation-vs-CONDITION-MATCHED-CONTROL. Controls are never pooled across
conditions, because IFNg stimulation and TIL co-culture shift the baseline
transcriptome enormously and pooling would manufacture effects that are really
just condition effects.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

__all__ = ["assign_pseudoreplicates", "make_pseudobulk", "group_sizes"]


def group_sizes(adata, cfg: dict[str, Any]) -> pd.DataFrame:
    """Cells per (perturbation x condition). The power table for the project."""
    s = cfg["schema"]["obs"]
    return (
        adata.obs.groupby([s["perturbation"], s["condition"]], observed=True)
        .size()
        .unstack(fill_value=0)
        .sort_index()
    )


def assign_pseudoreplicates(adata, cfg: dict[str, Any], key: str = "pseudorep") -> None:
    """Randomly split each perturbation x condition group into N pseudo-reps.

    CAVEAT, and state it plainly in the README: these are technical splits of
    one biological sample, not biological replicates. They give DESeq2 a
    within-group variance term to work with, but they systematically
    UNDERSTATE biological variance, so p-values are anti-conservative. The
    E-distance permutation test in nb03 is the honest effect-size check;
    treat DE ranks as a ranking, not as calibrated inference.
    """
    s = cfg["schema"]["obs"]
    n_rep = int(cfg["de"]["n_pseudoreplicates"])
    rng = np.random.default_rng(int(cfg["project"]["seed"]))

    labels = pd.Series(index=adata.obs_names, dtype=object)
    grouped = adata.obs.groupby([s["perturbation"], s["condition"]], observed=True)
    for (pert, cond), idx in grouped.groups.items():
        idx = np.asarray(idx)
        draw = rng.integers(0, n_rep, size=len(idx))
        labels.loc[idx] = [f"{pert}|{cond}|r{d}" for d in draw]
    adata.obs[key] = labels.astype("category")


def make_pseudobulk(adata, cfg: dict[str, Any], key: str = "pseudorep", layer: str | None = None):
    """Sum raw counts within each pseudo-replicate.

    Returns (counts_df [reps x genes], meta_df [reps x covariates]).
    Summing COUNTS (not normalised values) is what makes the result a valid
    input to a negative-binomial model such as DESeq2/pydeseq2.
    """
    import scipy.sparse as sp

    X = adata.layers[layer] if layer else adata.X
    groups = adata.obs[key].astype(str).values
    uniq = pd.unique(groups)
    index = {g: i for i, g in enumerate(uniq)}

    # One-hot design (reps x cells) @ (cells x genes) -> (reps x genes)
    rows = np.array([index[g] for g in groups])
    cols = np.arange(len(groups))
    M = sp.csr_matrix(
        (np.ones(len(groups)), (rows, cols)), shape=(len(uniq), adata.n_obs)
    )
    summed = M @ X
    summed = np.asarray(summed.todense()) if sp.issparse(summed) else np.asarray(summed)

    counts = pd.DataFrame(
        np.rint(summed).astype(int), index=uniq, columns=adata.var_names
    )

    parts = pd.Series(uniq).str.split("|", expand=True)
    meta = pd.DataFrame(
        {"perturbation": parts[0].values, "condition": parts[1].values, "rep": parts[2].values},
        index=uniq,
    )
    meta["n_cells"] = pd.Series(groups).value_counts().reindex(uniq).values

    keep = meta["n_cells"] >= int(cfg["de"]["min_cells_per_pseudoreplicate"])
    if (~keep).any():
        print(f"[filter] dropping {(~keep).sum()} pseudo-reps below min cell count")
    return counts.loc[keep], meta.loc[keep]
