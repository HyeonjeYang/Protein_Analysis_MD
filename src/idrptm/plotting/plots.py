"""Figure-generation placeholders."""

from __future__ import annotations

from pathlib import Path


def save_placeholder_plot(output: str | Path) -> Path:
    """Reserve a plotting API boundary for Stage 2."""

    path = Path(output)
    raise NotImplementedError(f"Plot generation will be implemented in Stage 2: {path}")
