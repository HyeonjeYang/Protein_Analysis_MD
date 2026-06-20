"""Residue parameter metadata placeholders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ResidueParameterSet:
    """Reference to a CALVADOS-compatible residue parameter table."""

    name: str
    path: Path | None = None
    supports_ptms: tuple[str, ...] = ()


def default_parameter_sets() -> list[ResidueParameterSet]:
    """Return known parameter-set placeholders.

    Stage 2 will wire these to concrete CALVADOS-compatible residue tables.
    """

    return [
        ResidueParameterSet(name="CALVADOS2", supports_ptms=()),
        ResidueParameterSet(name="pCALVADOS2", supports_ptms=("pSer", "pThr")),
    ]
