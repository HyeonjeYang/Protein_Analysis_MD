"""Public package namespace for protein_analysis_md.

The historical ``idrptm`` package remains available as a backward-compatible
namespace while the project migrates to ``protein_analysis_md``.
"""

from __future__ import annotations

import importlib
import sys

from idrptm import SUPPORTED_MVP_PTMS, __version__

_ALIASES = [
    "analysis",
    "analysis.cleavage",
    "analysis.compare",
    "analysis.contacts",
    "analysis.energy",
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
    "calvados_adapter",
    "cleavage",
    "design",
    "enzymes",
    "environment",
    "hpc",
    "ptm",
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
]

for _name in _ALIASES:
    try:
        sys.modules[f"{__name__}.{_name}"] = importlib.import_module(f"idrptm.{_name}")
    except ModuleNotFoundError:
        continue

__all__ = ["SUPPORTED_MVP_PTMS", "__version__"]
