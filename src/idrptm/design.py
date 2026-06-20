"""Variant design placeholders."""

from __future__ import annotations

from dataclasses import dataclass

from idrptm.ptm import PTMState
from idrptm.sequence import SequenceRecord


@dataclass(frozen=True)
class DesignedVariant:
    """A sequence/PTM design ready to be translated into backend input."""

    name: str
    sequence: SequenceRecord
    ptm_state: PTMState


def design_variants(sequence: SequenceRecord) -> list[DesignedVariant]:
    """Return the Stage 1 WT-only placeholder design."""

    wt_state = PTMState(name="WT")
    return [DesignedVariant(name=f"{sequence.name}_WT", sequence=sequence, ptm_state=wt_state)]
