"""Simple equilibration diagnostics for trajectory-derived time series."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def running_mean(values: np.ndarray, window: int) -> np.ndarray:
    """Compute a trailing running mean."""

    if window <= 1:
        return values.astype(float)
    series = pd.Series(values, dtype=float)
    return series.rolling(window=window, min_periods=1).mean().to_numpy()


def block_means(values: np.ndarray, n_blocks: int = 5) -> pd.DataFrame:
    """Compute block means for a one-dimensional time series."""

    if values.size == 0:
        return pd.DataFrame(columns=["block", "start_frame", "end_frame", "mean"])
    blocks = np.array_split(np.arange(values.size), min(n_blocks, values.size))
    rows = []
    for block_index, frame_indices in enumerate(blocks):
        rows.append(
            {
                "block": block_index,
                "start_frame": int(frame_indices[0]),
                "end_frame": int(frame_indices[-1]),
                "mean": float(np.mean(values[frame_indices])),
            }
        )
    return pd.DataFrame(rows)


def equilibration_diagnostics(
    *,
    rg: np.ndarray | None = None,
    potential_energy: np.ndarray | None = None,
    window: int | None = None,
) -> dict[str, object]:
    """Return cautious stability diagnostics from Rg and/or potential energy."""

    series = potential_energy if potential_energy is not None else rg
    if series is None or len(series) == 0:
        return {
            "available": False,
            "message": "No Rg or potential-energy time series available.",
        }
    values = np.asarray(series, dtype=float)
    run_window = window or max(1, len(values) // 10)
    means = running_mean(values, run_window)
    first = float(np.mean(values[: max(1, len(values) // 5)]))
    last = float(np.mean(values[-max(1, len(values) // 5) :]))
    scale = float(np.std(values)) or 1.0
    stable = abs(last - first) <= 2.0 * scale
    discard = int(0.2 * len(values)) if not stable else int(0.1 * len(values))
    return {
        "available": True,
        "metric": "potential_energy" if potential_energy is not None else "rg",
        "running_mean_window": run_window,
        "first_block_mean": first,
        "last_block_mean": last,
        "std": scale,
        "simple_stability_heuristic": "stable" if stable else "drifting",
        "recommended_discard_frames": discard,
        "caution": "This is a heuristic diagnostic, not proof of equilibrium.",
        "running_mean_last": float(means[-1]),
    }


def write_equilibration_outputs(analysis_dir: str | Path) -> dict[str, Path]:
    """Write diagnostics JSON and block means from existing analysis outputs."""

    root = Path(analysis_dir)
    rg_path = root / "timeseries_rg.parquet"
    energy_path = root / "energy.parquet"
    rg_values = None
    energy_values = None
    if rg_path.exists():
        rg_values = pd.read_parquet(rg_path)["rg"].to_numpy(dtype=float)
    if energy_path.exists():
        table = pd.read_parquet(energy_path)
        if "potential_energy_kj_mol" in table:
            energy_values = table["potential_energy_kj_mol"].to_numpy(dtype=float)
    diagnostics = equilibration_diagnostics(rg=rg_values, potential_energy=energy_values)
    json_path = root / "equilibration_diagnostics.json"
    json_path.write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    blocks_path = root / "equilibration_blocks.parquet"
    values = energy_values if energy_values is not None else rg_values
    if values is None:
        block_means(np.array([])).to_parquet(blocks_path, index=False)
    else:
        block_means(np.asarray(values, dtype=float)).to_parquet(blocks_path, index=False)
    return {"equilibration_diagnostics": json_path, "equilibration_blocks": blocks_path}
