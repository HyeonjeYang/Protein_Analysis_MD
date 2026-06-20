"""Sequence-feature analysis placeholders."""

from __future__ import annotations

from idrptm.sequence import SequenceRecord


def amino_acid_composition(sequence: SequenceRecord) -> dict[str, float]:
    """Return simple residue fractions for early sanity checks."""

    return {
        residue: sequence.sequence.count(residue) / sequence.length
        for residue in sorted(set(sequence.sequence))
    }
