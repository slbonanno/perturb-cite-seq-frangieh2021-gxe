"""Download, load and schema-check the Frangieh 2021 Perturb-CITE-seq data.

Named `data_io` rather than `io` to avoid shadowing the stdlib module.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from .config import load_config, paths

__all__ = [
    "download_file",
    "fetch_raw",
    "load_rna",
    "load_protein",
    "describe_schema",
    "check_schema",
    "align_modalities",
]


def download_file(url: str, dest: Path, overwrite: bool = False) -> Path:
    """Stream a URL to disk, skipping if already present."""
    dest = Path(dest)
    if dest.exists() and not overwrite:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size / 1e9:.2f} GB)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[get ] {dest.name} <- {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    print(f"[done] {dest.name} ({dest.stat().st_size / 1e9:.2f} GB)")
    return dest


def fetch_raw(cfg: dict[str, Any] | None = None, overwrite: bool = False) -> dict[str, Path]:
    """Download both modalities into data/raw/ per config."""
    cfg = cfg or load_config()
    P = paths(cfg)
    out = {}
    for key, spec in cfg["data"]["files"].items():
        out[key] = download_file(spec["url"], P.data_raw / spec["filename"], overwrite)
    return out


def _read(path: Path, backed: str | None = None):
    import anndata as ad

    return ad.read_h5ad(path, backed=backed)


def load_rna(cfg: dict[str, Any] | None = None, backed: str | None = None):
    cfg = cfg or load_config()
    P = paths(cfg)
    return _read(P.data_raw / cfg["data"]["files"]["rna"]["filename"], backed=backed)


def load_protein(cfg: dict[str, Any] | None = None):
    cfg = cfg or load_config()
    P = paths(cfg)
    return _read(P.data_raw / cfg["data"]["files"]["protein"]["filename"])


# -----------------------------------------------------------------------------
# Schema discovery -- run this FIRST, in nb01, before writing any analysis.
# -----------------------------------------------------------------------------

def describe_schema(adata, name: str = "adata", max_levels: int = 12) -> pd.DataFrame:
    """Print and return a summary of every .obs column.

    The config `schema:` block is a guess until this has been run. Use the
    output to correct config/config.yaml, then re-run everything downstream.
    """
    rows = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        nuniq = s.nunique(dropna=True)
        levels = (
            ", ".join(map(str, sorted(s.dropna().unique())[:max_levels]))
            if nuniq <= max_levels
            else f"<{nuniq} distinct values>"
        )
        rows.append(
            {"column": col, "dtype": str(s.dtype), "n_unique": nuniq, "levels": levels}
        )
    df = pd.DataFrame(rows)
    print(f"--- {name}: {adata.n_obs:,} cells x {adata.n_vars:,} features ---")
    print(f".obs columns ({len(adata.obs.columns)}):")
    print(df.to_string(index=False))
    print(f"\n.var columns: {list(adata.var.columns)}")
    print(f".obsm keys:   {list(adata.obsm.keys())}")
    print(f".layers keys: {list(adata.layers.keys())}")
    return df


def check_schema(adata, cfg: dict[str, Any] | None = None, strict: bool = True) -> bool:
    """Assert that the columns named in config actually exist.

    Call this at the top of every notebook after loading. Failing loudly here
    is much cheaper than a silent KeyError forty cells later.
    """
    cfg = cfg or load_config()
    want = cfg["schema"]["obs"]
    missing = {k: v for k, v in want.items() if v not in adata.obs.columns}
    if missing:
        msg = (
            "Config schema does not match the data.\n"
            f"  Missing .obs columns: {missing}\n"
            f"  Available: {list(adata.obs.columns)}\n"
            "Fix config/config.yaml -> schema.obs, then re-run."
        )
        if strict:
            raise KeyError(msg)
        print("[warn] " + msg)
        return False

    # Control label and condition levels
    pert_col = want["perturbation"]
    cond_col = want["condition"]
    ctrl = cfg["schema"]["control_label"]
    if ctrl not in set(adata.obs[pert_col].astype(str)):
        raise ValueError(
            f"control_label {ctrl!r} not found in obs[{pert_col!r}]. "
            f"Observed examples: {sorted(set(adata.obs[pert_col].astype(str)))[:10]}"
        )
    observed = sorted(set(adata.obs[cond_col].astype(str)))
    declared = cfg["schema"]["conditions"]["levels"]
    if set(observed) != set(declared):
        print(
            f"[warn] condition levels differ.\n"
            f"       config:   {declared}\n"
            f"       observed: {observed}"
        )
    print("[ok  ] schema check passed")
    return True


def align_modalities(rna, adt):
    """Intersect RNA and ADT on shared cell barcodes, in matching order.

    The two h5ads are distributed separately and are NOT guaranteed to contain
    the same cells in the same order. Always align before pairing modalities.
    """
    shared = rna.obs_names.intersection(adt.obs_names)
    print(
        f"RNA {rna.n_obs:,} | ADT {adt.n_obs:,} | shared {len(shared):,} "
        f"({len(shared) / max(rna.n_obs, 1):.1%} of RNA)"
    )
    return rna[shared].copy(), adt[shared].copy()
