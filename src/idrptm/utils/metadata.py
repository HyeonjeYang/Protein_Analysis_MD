"""Metadata JSON helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_metadata(path: str | Path, metadata: dict[str, Any]) -> Path:
    """Write metadata as pretty JSON."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return output
