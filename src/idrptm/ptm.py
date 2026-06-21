"""PTM registry and application for the MVP phosphorylation scope."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pandas as pd

SupportedPTM = Literal["pSer", "pThr"]
SUPPORTED_PTMS: tuple[SupportedPTM, ...] = ("pSer", "pThr")


@dataclass(frozen=True)
class PTMDefinition:
    """Definition of a supported residue-level PTM encoding."""

    name: str
    source_residue: str
    simulation_code: str
    description: str


@dataclass(frozen=True)
class PTMRequest:
    """A requested PTM using one-based biological residue numbering."""

    biological_position: int
    ptm: SupportedPTM
    expected_residue: str | None = None

    @property
    def zero_based_index(self) -> int:
        """Return the internal zero-based residue index."""

        return self.biological_position - 1


@dataclass(frozen=True)
class AppliedPTM:
    """A successfully applied PTM with both numbering conventions."""

    biological_position: int
    zero_based_index: int
    ptm: SupportedPTM
    source_residue: str
    simulation_code: str


@dataclass(frozen=True)
class PTMState:
    """A named PTM state for a sequence variant."""

    name: str
    sites: tuple[AppliedPTM, ...] = ()

    @property
    def is_wild_type(self) -> bool:
        """Return true when no PTM sites are present."""

        return not self.sites


def validate_supported_ptm(ptm: str) -> SupportedPTM:
    """Validate that a PTM is inside the MVP scope."""

    if ptm not in SUPPORTED_PTMS:
        supported = ", ".join(SUPPORTED_PTMS)
        raise ValueError(f"Unsupported PTM {ptm!r}. MVP supports: {supported}.")
    return ptm  # type: ignore[return-value]


class PTMRegistry:
    """Registry of supported PTM encodings."""

    def __init__(self, definitions: tuple[PTMDefinition, ...]) -> None:
        self._definitions = {definition.name: definition for definition in definitions}

    def get(self, ptm: str) -> PTMDefinition:
        """Return a PTM definition or raise a clear error."""

        supported = validate_supported_ptm(ptm)
        return self._definitions[supported]

    def apply(
        self,
        original_sequence: str,
        requests: tuple[PTMRequest, ...],
    ) -> tuple[str, tuple[AppliedPTM, ...]]:
        """Apply requested PTMs to a sequence.

        Inputs use one-based biological positions; validation and mutation use
        zero-based indices internally.
        """

        sequence_chars = list(original_sequence.upper())
        applied: list[AppliedPTM] = []
        seen_indices: set[int] = set()

        for request in sorted(requests, key=lambda item: item.zero_based_index):
            definition = self.get(request.ptm)
            index = request.zero_based_index
            if index < 0 or index >= len(sequence_chars):
                raise ValueError(
                    f"PTM site {request.biological_position} is outside sequence "
                    f"length {len(sequence_chars)}."
                )
            if index in seen_indices:
                raise ValueError(
                    f"Multiple PTMs requested for residue {request.biological_position}."
                )
            observed = sequence_chars[index]
            if request.expected_residue is not None and request.expected_residue != observed:
                raise ValueError(
                    f"Configured residue {request.expected_residue}{request.biological_position} "
                    f"does not match sequence residue {observed}{request.biological_position}."
                )
            if observed != definition.source_residue:
                raise ValueError(
                    f"{request.ptm} requires {definition.source_residue} at position "
                    f"{request.biological_position}, found {observed}."
                )
            sequence_chars[index] = definition.simulation_code
            seen_indices.add(index)
            applied.append(
                AppliedPTM(
                    biological_position=request.biological_position,
                    zero_based_index=index,
                    ptm=request.ptm,
                    source_residue=definition.source_residue,
                    simulation_code=definition.simulation_code,
                )
            )

        return "".join(sequence_chars), tuple(applied)


DEFAULT_PTM_REGISTRY = PTMRegistry(
    (
        PTMDefinition(
            name="pSer",
            source_residue="S",
            simulation_code="B",
            description="Phosphoserine encoded for CALVADOS-style simulation input.",
        ),
        PTMDefinition(
            name="pThr",
            source_residue="T",
            simulation_code="O",
            description="Phosphothreonine encoded for CALVADOS-style simulation input.",
        ),
    )
)


def apply_ptms(
    original_sequence: str,
    requests: tuple[PTMRequest, ...],
    registry: PTMRegistry = DEFAULT_PTM_REGISTRY,
) -> tuple[str, tuple[AppliedPTM, ...]]:
    """Apply pSer/pThr requests using the default MVP registry."""

    return registry.apply(original_sequence=original_sequence, requests=requests)


def parse_ptm_file(path: str | Path) -> tuple[PTMRequest, ...]:
    """Parse TXT/TSV/CSV PTM site files with line-numbered errors."""

    input_path = Path(path)
    requests: list[PTMRequest] = []
    lines = input_path.read_text(encoding="utf-8").splitlines()
    header: list[str] | None = None
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = _split_ptm_line(stripped)
        if header is None and {"site", "residue", "ptm"}.issubset({part.lower() for part in parts}):
            header = [part.lower() for part in parts]
            continue
        try:
            if header is not None:
                values = dict(zip(header, parts, strict=False))
                site = int(values["site"])
                residue = values["residue"]
                ptm = values["ptm"]
            else:
                site = int(parts[0])
                residue = parts[1]
                ptm = parts[2]
            requests.append(
                PTMRequest(
                    biological_position=site,
                    ptm=validate_supported_ptm(ptm),
                    expected_residue=residue,
                )
            )
        except Exception as exc:
            raise ValueError(
                f"Invalid PTM file {input_path} at line {line_number}: {line}"
            ) from exc
    return tuple(requests)


def build_ptm_table(sites: tuple[AppliedPTM, ...]) -> pd.DataFrame:
    """Build a PTM metadata table."""

    return pd.DataFrame(
        [
            {
                "site": site.biological_position,
                "zero_based_index": site.zero_based_index,
                "residue": site.source_residue,
                "ptm": site.ptm,
                "simulation_code": site.simulation_code,
            }
            for site in sites
        ],
        columns=["site", "zero_based_index", "residue", "ptm", "simulation_code"],
    )


def _split_ptm_line(line: str) -> list[str]:
    if "," in line:
        return [part.strip() for part in line.split(",")]
    return line.split()
