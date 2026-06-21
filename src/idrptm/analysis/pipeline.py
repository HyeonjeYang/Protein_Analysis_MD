"""Run the pure analysis core on loaded trajectory data and write artifacts."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances
from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.decomposition import run_decomposition_analysis
from idrptm.analysis.equilibration import write_equilibration_outputs
from idrptm.analysis.free_energy import run_free_energy_analysis
from idrptm.analysis.io import TrajectoryData, load_calvados_trajectory
from idrptm.analysis.lifetime import contact_lifetime
from idrptm.analysis.msd import com_msd
from idrptm.analysis.multichain import (
    chain_com_msd,
    cluster_size_timeseries,
    com_distance_timeseries,
    inter_chain_contact_map,
    inter_protein_contact_lifetime,
    intra_chain_contact_map,
    per_chain_ree,
    per_chain_rg,
)
from idrptm.analysis.ps import p_of_s
from idrptm.analysis.ree import ree_timeseries
from idrptm.analysis.rg import rg_timeseries
from idrptm.analysis.scaling import fit_flory_exponent, internal_distance_scaling
from idrptm.analysis.smoothing import (
    logspace_smooth_1d,
    rolling_smooth_1d,
    savgol_smooth_1d,
    smooth_contact_map,
)
from idrptm.schema import AnalysisConfig, WorkflowConfig, load_config
from idrptm.units import analysis_output_units, summary_units, write_units_metadata
from idrptm.visualization.smoothing_policy import validate_smoothing_request


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
    force: bool = False,
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
    return analyze_trajectory_data(
        data,
        output_dir=root,
        analysis_config=analysis_config,
        force=force,
    )


def analyze_trajectory_data(
    trajectory: TrajectoryData,
    *,
    output_dir: str | Path,
    analysis_config: AnalysisConfig | None = None,
    force: bool = False,
) -> AnalysisResult:
    """Run configured analyses for an already-loaded trajectory."""

    config = analysis_config or AnalysisConfig()
    enabled = {observable.lower() for observable in config.observables}
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    cache_key = _cache_key(trajectory, config)
    cached = _cached_result(output_path, cache_key)
    if cached is not None and not force:
        return cached
    positions = as_position_trajectory(trajectory.positions)
    frame_data = _frame_table(trajectory)
    outputs: dict[str, Path] = {}

    rg = rg_timeseries(positions)
    rg_table = frame_data.copy()
    rg_table["rg"] = rg
    _append_timeseries_smoothing(rg_table, config, key="rg", value_column="rg")
    outputs["timeseries_rg"] = _write_parquet(
        rg_table,
        output_path / "timeseries_rg.parquet",
        output_name="timeseries_rg",
    )

    ree = ree_timeseries(positions)
    ree_table = frame_data.copy()
    ree_table["ree"] = ree
    _append_timeseries_smoothing(ree_table, config, key="ree", value_column="ree")
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
    contact_map_smoothing = _enabled_smoothing(config, "contact_map")
    if contact_map_smoothing:
        smoothed_map = smooth_contact_map(
            contact_map,
            method=str(contact_map_smoothing.get("method", "none")),
            sigma=float(contact_map_smoothing.get("sigma", 1.0)),
            size=int(contact_map_smoothing.get("size", 3)),
            preserve_diagonal=bool(contact_map_smoothing.get("preserve_diagonal", True)),
        )
        smoothed_map_path = output_path / "contact_map_smoothed.npy"
        np.save(smoothed_map_path, smoothed_map)
        write_units_metadata(smoothed_map_path, analysis_output_units("contact_map"))
        outputs["contact_map_smoothed"] = smoothed_map_path

    ps_table = p_of_s(contact_map)
    _append_logspace_smoothing(ps_table, config, key="ps", value_column="p_contact")
    outputs["ps"] = _write_parquet(ps_table, output_path / "ps.parquet", output_name="ps")

    scaling = internal_distance_scaling(positions)
    _append_logspace_smoothing(
        scaling,
        config,
        key="rs",
        value_column="mean_distance_nm",
        output_column="mean_distance_nm_smooth",
    )
    outputs["scaling"] = _write_parquet(
        scaling,
        output_path / "scaling.parquet",
        output_name="scaling",
    )

    flory_fit = None
    if len(scaling) >= 2:
        if config.fit_to == "smoothed" and "mean_distance_nm_smooth" in scaling:
            fit_distances = scaling["mean_distance_nm_smooth"].to_numpy(dtype=float)
            flory_fit = fit_flory_exponent(
                s=scaling["s"].to_numpy(dtype=float),
                distances=fit_distances,
                min_s=config.fit_min_s,
                max_s=config.fit_max_s,
            )
        else:
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

    decomposition_outputs = run_decomposition_analysis(
        output_dir=output_path,
        positions=positions,
        frame_table=frame_data,
        contact_map=contact_map,
        analysis_config=config,
        rg=rg,
        ree=ree,
        sequence=_trajectory_sequence(trajectory),
    )
    outputs.update(decomposition_outputs)

    free_energy_outputs = run_free_energy_analysis(
        output_dir=output_path,
        frame_table=frame_data,
        analysis_config=config,
        rg=rg,
        ree=ree,
    )
    for path in free_energy_outputs.values():
        if path.suffix in {".csv", ".npy"}:
            write_units_metadata(path, analysis_output_units("free_energy"))
    outputs.update(free_energy_outputs)

    chain_ids = trajectory.residue_metadata()["chain_id"]
    if len(set(chain_ids.tolist())) > 1:
        outputs.update(
            _write_multichain_outputs(
                positions=positions,
                chain_ids=chain_ids,
                output_path=output_path,
                config=config,
                enabled=enabled,
            )
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
                smoothing=_summary_smoothing(config),
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["summary"] = summary_path
    write_equilibration_outputs(output_path)
    _write_cache_manifest(output_path, cache_key, outputs)
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
            "Install the declared protein_analysis_md dependencies."
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


def _write_multichain_outputs(
    *,
    positions: np.ndarray,
    chain_ids: np.ndarray,
    output_path: Path,
    config: AnalysisConfig,
    enabled: set[str],
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    outputs["per_chain_rg"] = _write_parquet(
        per_chain_rg(positions, chain_ids),
        output_path / "per_chain_rg.parquet",
        output_name="per_chain_rg",
    )
    outputs["per_chain_ree"] = _write_parquet(
        per_chain_ree(positions, chain_ids),
        output_path / "per_chain_ree.parquet",
        output_name="per_chain_ree",
    )
    intra_maps = intra_chain_contact_map(
        positions,
        chain_ids,
        cutoff=config.contact_cutoff_nm,
        min_sequence_separation=config.min_sequence_separation,
    )
    intra_path = output_path / "intra_chain_contact_map.npz"
    np.savez(intra_path, **intra_maps)
    write_units_metadata(intra_path, analysis_output_units("intra_chain_contact_map"))
    outputs["intra_chain_contact_map"] = intra_path

    inter_path = output_path / "inter_chain_contact_map.npy"
    np.save(
        inter_path,
        inter_chain_contact_map(
            positions,
            chain_ids,
            cutoff=config.contact_cutoff_nm,
        ),
    )
    write_units_metadata(inter_path, analysis_output_units("inter_chain_contact_map"))
    outputs["inter_chain_contact_map"] = inter_path

    outputs["com_distance"] = _write_parquet(
        com_distance_timeseries(positions, chain_ids),
        output_path / "com_distance.parquet",
        output_name="com_distance",
    )
    outputs["cluster_size"] = _write_parquet(
        cluster_size_timeseries(positions, chain_ids, cutoff=config.contact_cutoff_nm),
        output_path / "cluster_size.parquet",
        output_name="cluster_size",
    )
    if "msd" in enabled:
        outputs["chain_com_msd"] = _write_parquet(
            chain_com_msd(positions, chain_ids, max_lag=config.max_lag),
            output_path / "chain_com_msd.parquet",
            output_name="chain_com_msd",
        )
    if "lifetime" in enabled or "contact_lifetime" in enabled:
        outputs["inter_protein_contact_lifetime"] = _write_parquet(
            inter_protein_contact_lifetime(
                positions,
                chain_ids,
                cutoff=config.contact_cutoff_nm,
                max_lag=config.max_lag,
            ),
            output_path / "inter_protein_contact_lifetime.parquet",
            output_name="inter_protein_contact_lifetime",
        )
    return outputs


def _summary(
    *,
    trajectory: TrajectoryData,
    analysis_config: AnalysisConfig,
    rg: np.ndarray,
    ree: np.ndarray,
    contact_map: np.ndarray,
    flory_fit: object | None,
    outputs: dict[str, Path],
    smoothing: dict[str, object],
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
            "fit_to": analysis_config.fit_to,
            "decomposition": analysis_config.decomposition,
            "free_energy": getattr(analysis_config, "free_energy", {}),
        },
        "smoothing": smoothing,
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


def _enabled_smoothing(config: AnalysisConfig, key: str) -> dict[str, object]:
    smoothing = config.smoothing.get(key, {}) if isinstance(config.smoothing, dict) else {}
    if not isinstance(smoothing, dict) or not smoothing.get("enabled", False):
        return {}
    return validate_smoothing_request(key, smoothing)


def _append_logspace_smoothing(
    table: pd.DataFrame,
    config: AnalysisConfig,
    *,
    key: str,
    value_column: str,
    output_column: str | None = None,
) -> None:
    smoothing = _enabled_smoothing(config, key)
    if not smoothing:
        return
    method = str(smoothing.get("method", "logspace"))
    if method != "logspace":
        raise ValueError(f"{key} smoothing currently supports method='logspace'.")
    smoothed = logspace_smooth_1d(
        table["s"].to_numpy(dtype=float),
        table[value_column].to_numpy(dtype=float),
        window_log10=float(smoothing.get("window_log10", 0.2)),
        min_points=int(smoothing.get("min_points", 5)),
        robust=bool(smoothing.get("robust", True)),
    )
    target = output_column or f"{value_column}_smooth"
    table[target] = smoothed["y_smooth"].to_numpy(dtype=float)
    table["n_points_smooth_window"] = smoothed["n_points_window"].to_numpy(dtype=int)
    table["smoothing_method"] = smoothed["smoothing_method"].to_numpy()
    table["smoothing_window_log10"] = smoothed["smoothing_window_log10"].to_numpy(dtype=float)


def _append_timeseries_smoothing(
    table: pd.DataFrame,
    config: AnalysisConfig,
    *,
    key: str,
    value_column: str,
) -> None:
    smoothing = _enabled_smoothing(config, key)
    if not smoothing:
        return
    x_column = "time_ps" if "time_ps" in table else "frame"
    method = str(smoothing.get("method", "rolling"))
    if method == "rolling":
        smoothed = rolling_smooth_1d(
            table[x_column].to_numpy(dtype=float),
            table[value_column].to_numpy(dtype=float),
            window=int(smoothing.get("window", 25)),
            method=str(smoothing.get("rolling_method", "mean")),
            center=bool(smoothing.get("center", True)),
        )
    elif method == "savgol":
        smoothed = savgol_smooth_1d(
            table[x_column].to_numpy(dtype=float),
            table[value_column].to_numpy(dtype=float),
            window_length=int(smoothing.get("window_length", 25)),
            polyorder=int(smoothing.get("polyorder", 2)),
        )
    else:
        raise ValueError(f"{key} smoothing supports method='rolling' or 'savgol'.")
    table[f"{value_column}_nm_smooth"] = smoothed["y_smooth"].to_numpy(dtype=float)
    table[f"{value_column}_smoothing_method"] = smoothed["smoothing_method"].to_numpy()


def _summary_smoothing(config: AnalysisConfig) -> dict[str, object]:
    metadata: dict[str, object] = {
        "policy": "conservative",
        "raw_data_preserved": True,
        "quantitative_metrics_use_raw_by_default": config.fit_to == "raw",
        "contact_map_smoothing_default_enabled": False,
        "delta_contact_map_smoothing_default_enabled": False,
    }
    for key, value in (config.smoothing or {}).items():
        if isinstance(value, dict):
            metadata[key] = validate_smoothing_request(key, value)
    if "ps" in metadata:
        metadata["ps"] |= {
            "raw_column": "p_contact",
            "smoothed_column": "p_contact_smooth",
        }
    if "rs" in metadata:
        metadata["rs"] |= {
            "raw_column": "mean_distance_nm",
            "smoothed_column": "mean_distance_nm_smooth",
        }
    if "rg" in metadata:
        metadata["rg"] |= {"raw_column": "rg", "smoothed_column": "rg_nm_smooth"}
    if "ree" in metadata:
        metadata["ree"] |= {"raw_column": "ree", "smoothed_column": "ree_nm_smooth"}
    return metadata


def _trajectory_sequence(trajectory: TrajectoryData) -> str | None:
    """Return a one-letter sequence if future trajectory metadata provides one."""

    metadata = trajectory.residue_metadata()
    residue_names = metadata.get("residue_name")
    if residue_names is None:
        return None
    sequence = "".join(str(residue) for residue in residue_names)
    return sequence if len(sequence) == trajectory.n_residues else None


def _cache_key(trajectory: TrajectoryData, config: AnalysisConfig) -> dict[str, object]:
    return {
        "topology": _file_fingerprint(trajectory.topology_path),
        "trajectory": _file_fingerprint(trajectory.trajectory_path),
        "analysis_config": config.model_dump(mode="json"),
        "software": "protein_analysis_md",
        "units_version": 1,
    }


def _file_fingerprint(path: str | Path | None) -> dict[str, object]:
    if path is None:
        return {}
    file_path = Path(path)
    if not file_path.exists():
        return {
            "path": str(file_path),
            "exists": False,
        }
    stat = file_path.stat()
    return {
        "path": str(file_path),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _cached_result(output_path: Path, cache_key: dict[str, object]) -> AnalysisResult | None:
    cache_path = output_path / "cache_manifest.json"
    summary_path = output_path / "summary.json"
    if not cache_path.exists() or not summary_path.exists():
        return None
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    if payload.get("cache_key") != cache_key:
        return None
    outputs = {
        name: Path(path)
        for name, path in payload.get("outputs", {}).items()
        if isinstance(path, str) and os.path.exists(path)
    }
    outputs["summary"] = summary_path
    return AnalysisResult(output_dir=output_path, summary_json=summary_path, outputs=outputs)


def _write_cache_manifest(
    output_path: Path,
    cache_key: dict[str, object],
    outputs: dict[str, Path],
) -> None:
    cache_path = output_path / "cache_manifest.json"
    cache_path.write_text(
        json.dumps(
            {
                "cache_key": cache_key,
                "outputs": {name: str(path) for name, path in outputs.items()},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
