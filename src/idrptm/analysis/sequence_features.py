"""Sequence-feature helpers for quick design sanity checks."""

from __future__ import annotations

from idrptm.sequence import SequenceRecord


def amino_acid_composition(sequence: SequenceRecord | str) -> dict[str, float]:
    """Return simple residue fractions for early sanity checks."""

    sequence_text = sequence.sequence if isinstance(sequence, SequenceRecord) else sequence
    return {
        residue: sequence_text.count(residue) / len(sequence_text)
        for residue in sorted(set(sequence_text))
    }


def compute_sequence_features(sequence: SequenceRecord | str) -> dict[str, float | int]:
    """Compute lightweight sequence descriptors without requiring trajectory data."""

    sequence_text = sequence.sequence if isinstance(sequence, SequenceRecord) else sequence
    length = len(sequence_text)
    if length == 0:
        raise ValueError("Cannot compute sequence features for an empty sequence.")
    charged = set("DEKR")
    positive = set("KR")
    negative = set("DE")
    aromatic = set("FYW")
    composition = amino_acid_composition(sequence_text)
    return {
        "length": length,
        "fraction_charged": sum(sequence_text.count(residue) for residue in charged) / length,
        "fraction_positive": sum(sequence_text.count(residue) for residue in positive) / length,
        "fraction_negative": sum(sequence_text.count(residue) for residue in negative) / length,
        "fraction_aromatic": sum(sequence_text.count(residue) for residue in aromatic) / length,
        "fraction_proline": sequence_text.count("P") / length,
        "fraction_glycine": sequence_text.count("G") / length,
        **{f"fraction_{residue}": value for residue, value in composition.items()},
    }
