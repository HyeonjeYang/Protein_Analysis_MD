from __future__ import annotations

import csv
import json

import numpy as np
import pandas as pd
import pytest

from idrptm.analysis.compare import compare_project
from idrptm.analysis.energy import smooth_energy_timeseries
from idrptm.analysis.io import TrajectoryData
from idrptm.analysis.pipeline import analyze_trajectory_data
from idrptm.analysis.smoothing import (
    coarse_bin_curve,
    logspace_smooth_1d,
    rolling_smooth_1d,
    savgol_smooth_1d,
    smooth_contact_map,
)
from idrptm.plotting.plots import plot_lines
from idrptm.schema import AnalysisConfig


def test_logspace_smooth_1d_preserves_raw_values_and_length() -> None:
    x = np.arange(1, 8, dtype=float)
    y = np.array([1.0, 2.0, np.nan, 4.0, 8.0, 16.0, 32.0])

    result = logspace_smooth_1d(x, y, window_log10=0.5, min_points=1)

    np.testing.assert_allclose(result["x"], x)
    np.testing.assert_allclose(result["y_raw"], y, equal_nan=True)
    assert len(result["y_smooth"]) == len(x)
    assert result["y_smooth"].notna().any()


def test_logspace_smooth_1d_requires_positive_x() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        logspace_smooth_1d([0.0, 1.0], [1.0, 2.0])


def test_rolling_smooth_1d_works_on_simple_series() -> None:
    result = rolling_smooth_1d([0, 1, 2], [1.0, 3.0, 5.0], window=2, center=False)

    assert result["y_smooth"].tolist() == pytest.approx([1.0, 2.0, 4.0])
    assert result["smoothing_method"].tolist() == ["rolling_mean"] * 3


def test_savgol_smooth_1d_validates_odd_window_length() -> None:
    with pytest.raises(ValueError, match="odd"):
        savgol_smooth_1d([0, 1, 2, 3], [1.0, 2.0, 3.0, 4.0], window_length=4, polyorder=2)


def test_contact_map_and_coarse_binning_smoothing_helpers() -> None:
    matrix = np.array([[1.0, 0.2], [0.4, 1.0]])
    smoothed = smooth_contact_map(matrix, method="gaussian", sigma=0.5)
    binned = coarse_bin_curve([1, 2, 4, 8], [1.0, 2.0, 4.0, 8.0], bins="log", n_bins=2)

    assert smoothed.shape == matrix.shape
    np.testing.assert_allclose(smoothed, smoothed.T)
    np.testing.assert_allclose(np.diag(smoothed), np.diag(matrix))
    assert binned["n_points"].sum() == 4


def test_energy_smoothing_appends_smoothed_columns() -> None:
    table = pd.DataFrame(
        {
            "time_ns": [0.0, 1.0, 2.0],
            "potential_energy_kj_mol": [1.0, 3.0, 5.0],
            "temperature_K": [300.0, 302.0, 304.0],
        }
    )

    smoothed = smooth_energy_timeseries(
        table,
        {"enabled": True, "method": "rolling", "window": 2, "center": False},
    )

    assert smoothed["potential_energy_kj_mol_smooth"].tolist() == pytest.approx(
        [1.0, 2.0, 4.0]
    )
    assert smoothed["temperature_K_smooth"].tolist() == pytest.approx([300.0, 301.0, 303.0])


def test_analysis_pipeline_writes_smoothed_outputs_and_metadata(tmp_path) -> None:
    trajectory = TrajectoryData(
        positions=np.array(
            [
                [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]],
                [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0], [1.8, 0.0, 0.0], [3.5, 0.0, 0.0]],
                [[0.2, 0.0, 0.0], [0.9, 0.0, 0.0], [2.0, 0.0, 0.0], [3.6, 0.0, 0.0]],
            ]
        ),
        topology_path="top.pdb",
        trajectory_path="trajectory.dcd",
        time_ps=np.array([0.0, 10.0, 20.0]),
    )
    config = AnalysisConfig(
        observables=["rg", "ree", "contacts", "ps", "scaling"],
        smoothing={
            "ps": {"enabled": True, "method": "logspace", "window_log10": 0.6, "min_points": 1},
            "rs": {"enabled": True, "method": "logspace", "window_log10": 0.6, "min_points": 1},
            "rg": {"enabled": True, "method": "rolling", "window": 2},
            "ree": {"enabled": True, "method": "rolling", "window": 2},
            "contact_map": {"enabled": True, "method": "gaussian", "sigma": 0.5},
        },
    )

    result = analyze_trajectory_data(
        trajectory,
        output_dir=tmp_path / "analysis",
        analysis_config=config,
    )

    ps = pd.read_parquet(result.output_dir / "ps.parquet")
    scaling = pd.read_parquet(result.output_dir / "scaling.parquet")
    rg = pd.read_parquet(result.output_dir / "timeseries_rg.parquet")
    ree = pd.read_parquet(result.output_dir / "timeseries_ree.parquet")
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert "p_contact_smooth" in ps
    assert "mean_distance_nm_smooth" in scaling
    assert "rg_nm_smooth" in rg
    assert "ree_nm_smooth" in ree
    assert (result.output_dir / "contact_map_smoothed.npy").is_file()
    assert summary["smoothing"]["ps"]["raw_column"] == "p_contact"
    assert summary["smoothing"]["rs"]["smoothed_column"] == "mean_distance_nm_smooth"
    assert summary["outputs"]["contact_map_smoothed"].endswith("contact_map_smoothed.npy")


def test_smoothed_contact_map_is_not_used_for_delta_metrics_by_default(tmp_path) -> None:
    _write_minimal_project(tmp_path)

    comparison = compare_project(tmp_path)

    assert comparison.delta_maps["pSer2"][0, 1] == pytest.approx(0.2)


def test_plot_lines_labels_smoothed_trend_when_used() -> None:
    table = pd.DataFrame(
        {
            "condition": ["WT", "WT", "WT"],
            "s": [1, 2, 3],
            "p_mean": [0.3, 0.2, 0.1],
            "p_smooth_mean": [0.28, 0.21, 0.12],
        }
    )

    fig = plot_lines(
        table,
        "s",
        "p_mean",
        "P(s) (dimensionless)",
        "Contact probability P(s)",
        smooth_y="p_smooth_mean",
    )

    assert "smoothed trend" in fig.axes[0].get_title()
    assert fig.axes[0].get_ylabel() == "P(s) (dimensionless)"


def _write_minimal_project(project_dir) -> None:
    rows = [
        {
            "variant_id": "toy__WT",
            "condition": "WT",
            "replicate": "WT_1",
            "ptm_state": "WT",
            "ptm_sites_1based": "",
            "original_sequence": "AST",
            "metadata_path": "runs/toy__WT/metadata.yaml",
        },
        {
            "variant_id": "toy__pSer2",
            "condition": "pSer2",
            "replicate": "rep1",
            "ptm_state": "pSer2",
            "ptm_sites_1based": "pSer:S2->B",
            "original_sequence": "AST",
            "metadata_path": "runs/toy__pSer2/metadata.yaml",
        },
    ]
    with (project_dir / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    _write_minimal_analysis(
        project_dir / "runs" / "toy__WT" / "analysis",
        contact_map=np.array([[0.0, 0.1], [0.1, 0.0]]),
        smoothed_contact_map=np.array([[0.0, 0.9], [0.9, 0.0]]),
    )
    _write_minimal_analysis(
        project_dir / "runs" / "toy__pSer2" / "analysis",
        contact_map=np.array([[0.0, 0.3], [0.3, 0.0]]),
        smoothed_contact_map=np.array([[0.0, 0.1], [0.1, 0.0]]),
    )


def _write_minimal_analysis(
    analysis_dir,
    *,
    contact_map: np.ndarray,
    smoothed_contact_map: np.ndarray,
) -> None:
    analysis_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"frame": [0, 1], "rg": [1.0, 1.1]}).to_parquet(
        analysis_dir / "timeseries_rg.parquet",
        index=False,
    )
    pd.DataFrame({"frame": [0, 1], "ree": [2.0, 2.1]}).to_parquet(
        analysis_dir / "timeseries_ree.parquet",
        index=False,
    )
    np.save(analysis_dir / "contact_map.npy", contact_map)
    np.save(analysis_dir / "contact_map_smoothed.npy", smoothed_contact_map)
    pd.DataFrame({"s": [1], "p": [float(contact_map[0, 1])], "n_pairs": [1]}).to_parquet(
        analysis_dir / "ps.parquet",
        index=False,
    )
    pd.DataFrame({"s": [1, 2], "distance": [1.0, 1.5], "n_pairs": [1, 1]}).to_parquet(
        analysis_dir / "scaling.parquet",
        index=False,
    )
    (analysis_dir / "summary.json").write_text(
        json.dumps({"flory_fit": {"nu": 0.5}}) + "\n",
        encoding="utf-8",
    )
