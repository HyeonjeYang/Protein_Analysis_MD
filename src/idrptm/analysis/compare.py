"""WT-vs-PTM comparison placeholders."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ComparisonResult:
    """Minimal named comparison result."""

    reference: str
    variant: str
    metric: str
    delta: float | None = None


def compare_observable(reference: str, variant: str, metric: str) -> ComparisonResult:
    """Create a placeholder comparison result."""

    return ComparisonResult(reference=reference, variant=variant, metric=metric)
