"""PTM state placeholders for the MVP phosphorylation scope."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

SupportedPTM = Literal["pSer", "pThr"]
SUPPORTED_PTMS: tuple[SupportedPTM, ...] = ("pSer", "pThr")


@dataclass(frozen=True)
class PTMState:
    """A named PTM state for a sequence variant."""

    name: str
    sites: tuple[int, ...] = ()
    ptm: SupportedPTM | None = None

    @property
    def is_wild_type(self) -> bool:
        """Return true when no PTM sites are present."""

        return not self.sites


def validate_supported_ptm(ptm: str) -> SupportedPTM:
    """Validate that a PTM is inside the Stage 1 MVP scope."""

    if ptm not in SUPPORTED_PTMS:
        supported = ", ".join(SUPPORTED_PTMS)
        raise ValueError(f"Unsupported PTM {ptm!r}. MVP supports: {supported}.")
    return ptm  # type: ignore[return-value]
