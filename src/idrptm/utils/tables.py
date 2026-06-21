"""Small table I/O helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_table(table: pd.DataFrame, path: str | Path) -> Path:
    """Write CSV or parquet based on file suffix."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".parquet":
        table.to_parquet(output, index=False)
    else:
        table.to_csv(output, index=False)
    return output
