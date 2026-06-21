from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from idrptm.analysis.io import TrajectoryData
from idrptm.analysis.pipeline import analyze_trajectory_data
from idrptm.schema import AnalysisConfig
from idrptm.visualization.smoothing_policy import (
    conservative_smoothing_defaults,
    event_schedules_are_never_smoothed,
    validate_smoothing_request,
)


def _positions() -> np.ndarray:
    frames = []
    for frame in range(4):
        frames.append([[i * 0.4 + 0.02 * frame, 0.03 * frame, 0.0] for i in range(6)])
    return np.asarray(frames, dtype=float)


def test_ps_rs_smoothing_columns_and_summary_metadata(tmp_path) -> None:
    config = AnalysisConfig(
        smoothing={
            "ps": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 1.0,
                "min_points": 1,
                "robust": True,
            },
            "rs": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 1.0,
                "min_points": 1,
                "robust": True,
            },
        }
    )
    result = analyze_trajectory_data(
        TrajectoryData(
            positions=_positions(),
            topology_path="top.pdb",
            trajectory_path="trajectory.dcd",
            time_ps=np.array([0.0, 10.0, 20.0, 30.0]),
        ),
        output_dir=tmp_path,
        analysis_config=config,
    )

    ps = pd.read_parquet(result.outputs["ps"])
    scaling = pd.read_parquet(result.outputs["scaling"])
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert "p_contact_smooth" in ps
    assert "mean_distance_nm_smooth" in scaling
    assert summary["smoothing"]["raw_data_preserved"] is True
    assert summary["smoothing"]["quantitative_metrics_use_raw_by_default"] is True


def test_summary_metrics_use_raw_values_by_default(tmp_path) -> None:
    result = analyze_trajectory_data(
        TrajectoryData(
            positions=_positions(),
            topology_path="top.pdb",
            trajectory_path="trajectory.dcd",
        ),
        output_dir=tmp_path,
        analysis_config=AnalysisConfig(
            smoothing={"rg": {"enabled": True, "method": "rolling", "window": 2}}
        ),
    )
    rg = pd.read_parquet(result.outputs["timeseries_rg"])
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert summary["rg_mean"] == pytest.approx(float(rg["rg"].mean()))


def test_contact_map_and_delta_smoothing_disabled_by_default() -> None:
    defaults = conservative_smoothing_defaults()

    assert defaults["contact_map"]["enabled"] is False
    assert defaults["delta_contact_map"]["enabled"] is False


def test_event_schedules_are_never_smoothed() -> None:
    assert event_schedules_are_never_smoothed()
    with pytest.raises(ValueError, match="not allowed"):
        validate_smoothing_request("event_schedule", {"enabled": True})


def test_smoothing_metadata_policy_for_contact_map_warns() -> None:
    metadata = validate_smoothing_request("contact_map", {"enabled": True, "method": "gaussian"})

    assert metadata["visualization_only"] is True
    assert metadata["warning_required"] is True
