"""Run the pure analysis core on loaded trajectory data and write artifacts."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances
from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.io import TrajectoryData, load_calvados_trajectory
from idrptm.analysis.lifetime import contact_lifetime
from idrptm.analysis.msd import com_msd
from idrptm.analysis.ps import p_of_s
from idrptm.analysis.ree import ree_timeseries
from idrptm.analysis.rg import rg_timeseries
from idrptm.analysis.scaling import fit_flory_exponent, internal_distance_scaling
from idrptm.schema import AnalysisConfig, WorkflowConfig, load_config
from idrptm.units import analysis_output_units, summary_units, write_units_metadata


@dataclass(frozen=True)
class AnalysisResult:
    """Paths written by the analyze workflow."""

    output_dir: Path
    summary_json: Path
    outputs: dict[str, Path]


def analyze_run_directory(
    run_dir: str | Path,
    *,
    config: WorkflowConfig | None = None,
    config_path: str | Path | None = None,
    topology: str | Path | None = None,
    trajectory: str | Path | None = None,
    trajectory_reader: Literal["mdtraj", "mdanalysis"] = "mdtraj",
    output_dir: str | Path | None = None,
) -> AnalysisResult:
    """Load a CALVADOS run directory and write analysis outputs."""

    workflow = config
    if workflow is None and config_path is not None:
        workflow = load_config(config_path)
    analysis_config = workflow.analysis if workflow is not None else AnalysisConfig()
    data = load_calvados_trajectory(
        run_dir,
        topology=topology,
        trajectory=trajectory,
        engine=trajectory_reader,
    )
    root = Path(output_dir) if output_dir is not None else Path(run_dir) / "analysis"
    return analyze_trajectory_data(data, output_dir=root, analysis_config=analysis_config)


def analyze_trajectory_data(
    trajectory: TrajectoryData,
    *,
    output_dir: str | Path,
    analysis_config: AnalysisConfig | None = None,
) -> AnalysisResult:
    """Run configured analyses for an already-loaded trajectory."""

    config = analysis_config or AnalysisConfig()
    enabled = {observable.lower() for observable in config.observables}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    positions = as_position_trajectory(trajectory.positions)
    frame_data = _frame_table(trajectory)
    outputs: dict[str, Path] = {}

    rg = rg_timeseries(positions)
    rg_table = frame_data.copy()
    rg_table["rg"] = rg
    outputs["timeseries_rg"] = _write_parquet(
        rg_table,
        output_path / "timeseries_rg.parquet",
        output_name="timeseries_rg",
    )

    ree = ree_timeseries(positions)
    ree_table = frame_data.copy()
    ree_table["ree"] = ree
    outputs["timeseries_ree"] = _write_parquet(
        ree_table,
        output_path / "timeseries_ree.parquet",
        output_name="timeseries_ree",
    )

    contact_map = contact_map_from_positions(
        positions,
        cutoff=config.contact_cutoff_nm,
        min_sequence_separation=config.min_sequence_separation,
    )
    contact_map_path = output_path / "contact_map.npy"
    np.save(contact_map_path, contact_map)
    write_units_metadata(contact_map_path, analysis_output_units("contact_map"))
    outputs["contact_map"] = contact_map_path

    ps_table = p_of_s(contact_map)
    outputs["ps"] = _write_parquet(ps_table, output_path / "ps.parquet", output_name="ps")

    scaling = internal_distance_scaling(positions)
    outputs["scaling"] = _write_parquet(
        scaling,
        output_path / "scaling.parquet",
        output_name="scaling",
    )

    flory_fit = None
    if len(scaling) >= 2:
        flory_fit = fit_flory_exponent(
            scaling,
            min_s=config.fit_min_s,
            max_s=config.fit_max_s,
        )

    if "msd" in enabled:
        outputs["msd"] = _write_parquet(
            com_msd(positions, max_lag=config.max_lag),
            output_path / "msd.parquet",
            output_name="msd",
        )

    if "lifetime" in enabled or "contact_lifetime" in enabled:
        contact_states = _contact_states(
            positions,
            cutoff=config.contact_cutoff_nm,
            min_sequence_separation=config.min_sequence_separation,
        )
        outputs["contact_lifetime"] = _write_parquet(
            contact_lifetime(contact_states, max_lag=config.max_lag),
            output_path / "contact_lifetime.parquet",
            output_name="contact_lifetime",
        )

    summary_path = output_path / "summary.json"
    summary_path.write_text(
        json.dumps(
            _summary(
                trajectory=trajectory,
                analysis_config=config,
                rg=rg,
                ree=ree,
                contact_map=contact_map,
                flory_fit=flory_fit,
                outputs=outputs,
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["summary"] = summary_path
    return AnalysisResult(output_dir=output_path, summary_json=summary_path, outputs=outputs)


def _frame_table(trajectory: TrajectoryData) -> pd.DataFrame:
    table = pd.DataFrame({"frame": np.arange(trajectory.n_frames, dtype=int)})
    if trajectory.time_ps is not None and len(trajectory.time_ps) == trajectory.n_frames:
        table["time_ps"] = trajectory.time_ps
    return table


def _write_parquet(table: pd.DataFrame, path: Path, *, output_name: str) -> Path:
    units = analysis_output_units(output_name)
    table.attrs["units"] = units
    try:
        table.to_parquet(path, index=False)
    except ImportError as exc:
        raise ImportError(
            "Writing parquet outputs requires pyarrow or fastparquet. "
            "Install the declared idr-ptm-md dependencies."
        ) from exc
    write_units_metadata(path, units)
    return path


def _contact_states(
    positions: np.ndarray,
    *,
    cutoff: float,
    min_sequence_separation: int,
) -> np.ndarray:
    n_residues = positions.shape[1]
    pair_i, pair_j = np.triu_indices(n_residues, k=min_sequence_separation)
    states = np.zeros((positions.shape[0], len(pair_i)), dtype=bool)
    for frame_index, frame in enumerate(positions):
        distances = pairwise_distances(frame)
        states[frame_index] = distances[pair_i, pair_j] <= cutoff
    return states


def _summary(
    *,
    trajectory: TrajectoryData,
    analysis_config: AnalysisConfig,
    rg: np.ndarray,
    ree: np.ndarray,
    contact_map: np.ndarray,
    flory_fit: object | None,
    outputs: dict[str, Path],
) -> dict[str, object]:
    contact_upper = contact_map[np.triu_indices_from(contact_map, k=1)]
    summary: dict[str, object] = {
        "topology_path": str(trajectory.topology_path),
        "trajectory_path": str(trajectory.trajectory_path),
        "length_unit": trajectory.length_unit,
        "input_position_unit": trajectory.input_position_unit,
        "canonical_position_unit": trajectory.canonical_position_unit,
        "n_frames": trajectory.n_frames,
        "n_residues": trajectory.n_residues,
        "units": summary_units(
            input_position_unit=trajectory.input_position_unit,
            canonical_position_unit=trajectory.canonical_position_unit,
        ),
        "analysis": {
            "observables": analysis_config.observables,
            "contact_cutoff_nm": analysis_config.contact_cutoff_nm,
            "min_sequence_separation": analysis_config.min_sequence_separation,
            "max_lag": analysis_config.max_lag,
        },
        "rg_mean": float(np.mean(rg)),
        "rg_std": float(np.std(rg)),
        "ree_mean": float(np.mean(ree)),
        "ree_std": float(np.std(ree)),
        "contact_probability_mean": float(contact_upper.mean()) if contact_upper.size else 0.0,
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    if flory_fit is not None:
        summary["flory_fit"] = {
            "nu": flory_fit.nu,
            "prefactor": flory_fit.prefactor,
            "r_value": flory_fit.r_value,
            "p_value": flory_fit.p_value,
            "stderr": flory_fit.stderr,
            "n_points": flory_fit.n_points,
        }
    return summary
