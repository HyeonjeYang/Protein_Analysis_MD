"""Sequence ingestion placeholders for Stage 1."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class SequenceRecord:
    """A single protein or IDR sequence record."""

    name: str
    sequence: str
    description: str = ""

    def __post_init__(self) -> None:
        sequence = self.sequence.replace(" ", "").replace("\n", "").upper()
        invalid = sorted(set(sequence) - CANONICAL_AMINO_ACIDS)
        if invalid:
            invalid_text = ", ".join(invalid)
            raise ValueError(f"Unsupported amino-acid code(s): {invalid_text}")
        if not sequence:
            raise ValueError("Sequence must not be empty.")
        object.__setattr__(self, "sequence", sequence)

    @property
    def length(self) -> int:
        """Sequence length in residues."""

        return len(self.sequence)


def read_fasta(path: str | Path) -> list[SequenceRecord]:
    """Read FASTA records.

    Stage 1 keeps this intentionally small; richer validation and sequence
    metadata handling belongs in Stage 2.
    """

    fasta_path = Path(path)
    records: list[SequenceRecord] = []
    name: str | None = None
    description = ""
    chunks: list[str] = []

    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if name is not None:
                records.append(
                    SequenceRecord(
                        name=name,
                        sequence="".join(chunks),
                        description=description,
                    )
                )
            header = stripped[1:]
            name, _, description = header.partition(" ")
            chunks = []
        else:
            chunks.append(stripped)

    if name is not None:
        records.append(
            SequenceRecord(name=name, sequence="".join(chunks), description=description)
        )
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}.")
    return records
