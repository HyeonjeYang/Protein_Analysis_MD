"""PTM registry and application for the MVP phosphorylation scope."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

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
