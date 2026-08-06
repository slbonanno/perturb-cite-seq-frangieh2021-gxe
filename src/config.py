"""Configuration loading and path resolution.

Every notebook starts with::

    from src.config import load_config, paths
    cfg = load_config()

Nothing downstream should hard-code a path, a threshold, or an .obs column
name. If you find yourself typing a magic number into a notebook, it belongs
in config/config.yaml instead.
"""

from __future__ import annotations

import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

__all__ = ["repo_root", "load_config", "load_panels", "paths", "set_seed", "Paths"]


def repo_root(start: Path | str | None = None) -> Path:
    """Walk upward until we find the directory containing config/config.yaml.

    Makes notebooks position-independent -- they work whether Jupyter was
    launched from the repo root or from notebooks/.
    """
    here = Path(start or Path.cwd()).resolve()
    for candidate in [here, *here.parents]:
        if (candidate / "config" / "config.yaml").is_file():
            return candidate
    raise FileNotFoundError(
        "Could not locate repo root (no config/config.yaml found walking up "
        f"from {here})."
    )


def load_config(path: Path | str | None = None) -> dict[str, Any]:
    """Load config/config.yaml as a plain dict."""
    root = repo_root()
    cfg_path = Path(path) if path else root / "config" / "config.yaml"
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["_root"] = str(root)
    return cfg


def load_panels(path: Path | str | None = None) -> dict[str, Any]:
    """Load config/panels.yaml as a plain dict."""
    root = repo_root()
    panel_path = Path(path) if path else root / "config" / "panels.yaml"
    with open(panel_path) as fh:
        return yaml.safe_load(fh)


@dataclass(frozen=True)
class Paths:
    """Resolved absolute paths for every project directory."""

    root: Path
    data_raw: Path
    data_interim: Path
    data_processed: Path
    figures: Path
    tables: Path

    def mkdirs(self) -> "Paths":
        for p in (
            self.data_raw,
            self.data_interim,
            self.data_processed,
            self.figures,
            self.tables,
        ):
            p.mkdir(parents=True, exist_ok=True)
        return self


def paths(cfg: dict[str, Any] | None = None) -> Paths:
    """Resolve the `project.paths` block into absolute paths."""
    cfg = cfg or load_config()
    root = Path(cfg["_root"])
    p = cfg["project"]["paths"]
    return Paths(
        root=root,
        data_raw=root / p["data_raw"],
        data_interim=root / p["data_interim"],
        data_processed=root / p["data_processed"],
        figures=root / p["figures"],
        tables=root / p["tables"],
    ).mkdirs()


def set_seed(cfg: dict[str, Any] | None = None) -> int:
    """Seed python, numpy and PYTHONHASHSEED from config."""
    cfg = cfg or load_config()
    seed = int(cfg["project"]["seed"])
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass
    return seed
