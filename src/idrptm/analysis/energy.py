"""Energy/state-log parsing for CALVADOS/OpenMM-style outputs."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pandas as pd

from idrptm.analysis.smoothing import rolling_smooth_1d, savgol_smooth_1d
from idrptm.units import analysis_output_units, write_units_metadata

ENERGY_SMOOTH_COLUMNS = (
    "potential_energy_kj_mol",
    "kinetic_energy_kj_mol",
    "total_energy_kj_mol",
    "temperature_K",
)


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


def smooth_energy_timeseries(
    table: pd.DataFrame,
    smoothing: dict[str, object] | None,
) -> pd.DataFrame:
    """Append optional smoothed energy and temperature columns."""

    if not smoothing or not smoothing.get("enabled", False):
        return table.copy()
    result = table.copy()
    method = str(smoothing.get("method", "rolling"))
    x_column = "time_ns" if "time_ns" in result else "step"
    for column in ENERGY_SMOOTH_COLUMNS:
        if column not in result:
            continue
        if method == "rolling":
            smoothed = rolling_smooth_1d(
                result[x_column],
                result[column],
                window=int(smoothing.get("window", 25)),
                method=str(smoothing.get("rolling_method", "mean")),
                center=bool(smoothing.get("center", True)),
            )
        elif method == "savgol":
            smoothed = savgol_smooth_1d(
                result[x_column],
                result[column],
                window_length=int(smoothing.get("window_length", 25)),
                polyorder=int(smoothing.get("polyorder", 2)),
            )
        else:
            raise ValueError("energy smoothing supports method='rolling' or 'savgol'.")
        result[f"{column}_smooth"] = smoothed["y_smooth"].to_numpy(dtype=float)
        result["energy_smoothing_method"] = smoothed["smoothing_method"].to_numpy()
    return result


def write_energy_outputs(
    table: pd.DataFrame,
    output_dir: str | Path,
    *,
    smoothing: dict[str, object] | None = None,
) -> tuple[Path, Path]:
    """Write ``energy.parquet`` and ``energy_summary.json``."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    parquet = root / "energy.parquet"
    summary = root / "energy_summary.json"
    output = smooth_energy_timeseries(table, smoothing)
    output.attrs["units"] = analysis_output_units("energy")
    output.to_parquet(parquet, index=False)
    write_units_metadata(parquet, analysis_output_units("energy"))
    summary.write_text(
        json.dumps(
            {
                "n_rows": int(len(output)),
                "units": {
                    "time_ns": "ns",
                    "energy": "kJ/mol",
                    "temperature": "K",
                },
                "smoothing": dict(smoothing or {}),
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
