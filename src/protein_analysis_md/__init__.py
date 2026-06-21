"""Public package namespace for protein_analysis_md.

The historical ``idrptm`` package remains available as a backward-compatible
namespace while the project migrates to ``protein_analysis_md``.
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))

from idrptm import SUPPORTED_MVP_PTMS, __version__

_ALIASES = [
    "analysis",
    "analysis.cleavage",
    "analysis.compare",
    "analysis.contacts",
    "analysis.decomposition",
    "analysis.energy",
    "analysis.equilibration",
    "analysis.free_energy",
    "analysis.io",
    "analysis.lifetime",
    "analysis.msd",
    "analysis.multichain",
    "analysis.phase",
    "analysis.pipeline",
    "analysis.ps",
    "analysis.ree",
    "analysis.rg",
    "analysis.scaling",
    "analysis.sequence_features",
    "analysis.smoothing",
    "bundle",
    "calvados_adapter",
    "cleavage",
    "configuration",
    "dashboard",
    "design",
    "enzymes",
    "environment",
    "environment_check",
    "finalize",
    "hpc",
    "local_pipeline",
    "presets",
    "project",
    "provenance",
    "ptm",
    "recipe",
    "registry",
    "repo_check",
    "residue_params",
    "runner",
    "schema",
    "sequence",
    "sequence_features",
    "storage",
    "uniprot",
    "units",
    "utils",
    "utils.logging",
    "utils.metadata",
    "utils.paths",
    "utils.tables",
    "visualization",
    "visualization.cleavage",
    "visualization.decomposition",
    "visualization.free_energy",
    "visualization.heatmaps",
    "visualization.phase",
    "visualization.pymol",
    "visualization.ptm",
    "visualization.sequence_tracks",
    "visualization.single_chain",
    "visualization.smoothing_policy",
    "watch",
]

for _name in _ALIASES:
    try:
        sys.modules[f"{__name__}.{_name}"] = importlib.import_module(f"idrptm.{_name}")
    except ModuleNotFoundError:
        continue

__all__ = ["SUPPORTED_MVP_PTMS", "__version__"]
