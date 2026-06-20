"""Sequence ingestion for raw and FASTA protein inputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
FASTA_WRAP = 80


def sanitize_identifier(name: str) -> str:
    """Return a stable filesystem-friendly identifier."""

    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    clean = clean.strip("._-")
    if not clean:
        raise ValueError("Identifier cannot be empty.")
    return clean


def normalize_raw_sequence(sequence: str) -> str:
    """Normalize a raw sequence string while preserving biological residues."""

    return re.sub(r"\s+", "", sequence).upper()


@dataclass(frozen=True)
class SequenceRecord:
    """A single protein or IDR sequence record."""

    name: str
    sequence: str
    description: str = ""

    def __post_init__(self) -> None:
        name = sanitize_identifier(self.name)
        sequence = normalize_raw_sequence(self.sequence)
        invalid = sorted(set(sequence) - CANONICAL_AMINO_ACIDS)
        if invalid:
            invalid_text = ", ".join(invalid)
            raise ValueError(f"Unsupported amino-acid code(s): {invalid_text}")
        if not sequence:
            raise ValueError("Sequence must not be empty.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "sequence", sequence)

    @property
    def length(self) -> int:
        """Sequence length in residues."""

        return len(self.sequence)


def parse_raw_sequence(name: str, sequence: str, description: str = "") -> SequenceRecord:
    """Parse an inline raw amino-acid sequence."""

    return SequenceRecord(name=name, sequence=sequence, description=description)


def read_fasta(path: str | Path) -> list[SequenceRecord]:
    """Read FASTA records.

    The parser accepts one or more records and keeps residue numbering aligned
    with the original sequence by rejecting gap characters and unsupported codes.
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


def load_single_sequence(
    *,
    name: str,
    sequence: str | None = None,
    fasta: str | Path | None = None,
) -> SequenceRecord:
    """Load exactly one sequence from either raw text or a FASTA file."""

    if bool(sequence) == bool(fasta):
        raise ValueError("Provide exactly one sequence source: sequence or fasta.")
    if sequence is not None:
        return parse_raw_sequence(name=name, sequence=sequence)

    assert fasta is not None
    records = read_fasta(fasta)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one FASTA record in {fasta}; found {len(records)}.")
    record = records[0]
    return SequenceRecord(
        name=name or record.name,
        sequence=record.sequence,
        description=record.description,
    )


def format_fasta(record: SequenceRecord, sequence: str | None = None) -> str:
    """Format a record as FASTA text."""

    output_sequence = record.sequence if sequence is None else normalize_raw_sequence(sequence)
    header = record.name
    if record.description:
        header = f"{header} {record.description}"
    lines = [f">{header}"]
    lines.extend(
        output_sequence[index : index + FASTA_WRAP]
        for index in range(0, len(output_sequence), FASTA_WRAP)
    )
    return "\n".join(lines) + "\n"


def write_fasta(record: SequenceRecord, path: str | Path, sequence: str | None = None) -> Path:
    """Write a single FASTA record and return the output path."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(format_fasta(record, sequence=sequence), encoding="utf-8")
    return output_path
