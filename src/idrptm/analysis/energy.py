"""Energy/state-log parsing for CALVADOS/OpenMM-style outputs."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pandas as pd


def parse_energy_log(path: str | Path, *, phase: str = "production") -> pd.DataFrame:
    """Parse a CSV/TSV energy log into canonical columns."""

    log_path = Path(path)
    if not log_path.is_file():
        warnings.warn(f"Missing energy log: {log_path}", stacklevel=2)
        return pd.DataFrame(
            columns=[
                "phase",
                "time_ns",
                "step",
                "potential_energy_kj_mol",
                "kinetic_energy_kj_mol",
                "total_energy_kj_mol",
                "temperature_K",
            ]
        )
    table = pd.read_csv(log_path, sep=None, engine="python", comment="#")
    table = table.rename(columns={column: _normalize_column(column) for column in table.columns})
    if "phase" not in table:
        table["phase"] = phase
    return table


def write_energy_outputs(table: pd.DataFrame, output_dir: str | Path) -> tuple[Path, Path]:
    """Write ``energy.parquet`` and ``energy_summary.json``."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    parquet = root / "energy.parquet"
    summary = root / "energy_summary.json"
    table.to_parquet(parquet, index=False)
    summary.write_text(
        json.dumps(
            {
                "n_rows": int(len(table)),
                "units": {
                    "time_ns": "ns",
                    "energy": "kJ/mol",
                    "temperature": "K",
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return parquet, summary


def _normalize_column(column: str) -> str:
    key = column.strip().lower().replace(" ", "_").replace("(", "").replace(")", "")
    aliases = {
        "time_ps": "time_ps",
        "time_ns": "time_ns",
        "step": "step",
        "potential_energy_kj/mol": "potential_energy_kj_mol",
        "potential_energy_kj_mol": "potential_energy_kj_mol",
        "kinetic_energy_kj/mol": "kinetic_energy_kj_mol",
        "kinetic_energy_kj_mol": "kinetic_energy_kj_mol",
        "total_energy_kj/mol": "total_energy_kj_mol",
        "total_energy_kj_mol": "total_energy_kj_mol",
        "temperature_k": "temperature_K",
    }
    return aliases.get(key, key)
