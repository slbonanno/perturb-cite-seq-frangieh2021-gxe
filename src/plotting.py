"""Shared plotting helpers. All styling comes from config -> plotting."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

__all__ = ["apply_style", "condition_palette", "savefig"]


def apply_style(cfg: dict[str, Any]) -> None:
    p = cfg["plotting"]
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": p["dpi"],
            "figure.figsize": tuple(p["figsize_default"]),
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "font.size": 10,
        }
    )


def condition_palette(cfg: dict[str, Any]) -> dict[str, str]:
    return dict(cfg["plotting"]["condition_colors"])


def savefig(fig, name: str, cfg: dict[str, Any], subdir: str = "") -> Path:
    """Write a figure to results/figures/<subdir>/<name>.<format>."""
    from .config import paths

    out = paths(cfg).figures / subdir
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.{cfg['plotting']['format']}"
    fig.savefig(path)
    print(f"[fig ] {path.relative_to(paths(cfg).root)}")
    return path
