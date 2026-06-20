from __future__ import annotations

import json
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from idrptm.analysis.io import TrajectoryData, TrajectoryFileError, load_calvados_trajectory
from idrptm.analysis.pipeline import analyze_trajectory_data
from idrptm.schema import AnalysisConfig


def _trajectory_data() -> TrajectoryData:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.75, 0.0, 0.0]],
            [[1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]],
        ]
    )
    return TrajectoryData(
        positions=positions,
        topology_path="top.pdb",
        trajectory_path="trajectory.dcd",
        time_ps=np.array([0.0, 10.0, 20.0]),
    )


def test_load_calvados_trajectory_uses_mdtraj_with_expected_files(tmp_path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    top = run_dir / "top.pdb"
    dcd = run_dir / "trajectory.dcd"
    top.write_text("fake top\n", encoding="utf-8")
    dcd.write_text("fake dcd\n", encoding="utf-8")
    calls = {}

    def fake_load(trajectory_path: str, top: str):
        calls["trajectory"] = trajectory_path
        calls["topology"] = top
        return SimpleNamespace(
            xyz=np.zeros((2, 3, 3), dtype=float),
            time=np.array([0.0, 5.0]),
        )

    monkeypatch.setitem(sys.modules, "mdtraj", SimpleNamespace(load=fake_load))

    data = load_calvados_trajectory(run_dir)

    assert calls == {"trajectory": str(dcd), "topology": str(top)}
    assert data.positions.shape == (2, 3, 3)
    assert data.topology_path == top
    assert data.trajectory_path == dcd
    assert data.time_ps.tolist() == [0.0, 5.0]


def test_load_calvados_trajectory_reports_missing_files(tmp_path) -> None:
    with pytest.raises(TrajectoryFileError, match="Missing topology PDB"):
        load_calvados_trajectory(tmp_path)

    (tmp_path / "top.pdb").write_text("fake top\n", encoding="utf-8")
    with pytest.raises(TrajectoryFileError, match="Missing trajectory DCD"):
        load_calvados_trajectory(tmp_path)


def test_analyze_trajectory_data_writes_requested_outputs(tmp_path) -> None:
    config = AnalysisConfig(
        observables=["rg", "ree", "contacts", "ps", "scaling", "msd", "lifetime"],
        contact_cutoff_nm=1.0,
        max_lag=1,
    )

    result = analyze_trajectory_data(
        _trajectory_data(),
        output_dir=tmp_path / "analysis",
        analysis_config=config,
    )

    expected = {
        "timeseries_rg.parquet",
        "timeseries_ree.parquet",
        "contact_map.npy",
        "ps.parquet",
        "scaling.parquet",
        "msd.parquet",
        "contact_lifetime.parquet",
        "summary.json",
    }
    assert {path.name for path in result.outputs.values()} == expected

    rg = pd.read_parquet(result.output_dir / "timeseries_rg.parquet")
    ree = pd.read_parquet(result.output_dir / "timeseries_ree.parquet")
    ps = pd.read_parquet(result.output_dir / "ps.parquet")
    scaling = pd.read_parquet(result.output_dir / "scaling.parquet")
    msd = pd.read_parquet(result.output_dir / "msd.parquet")
    lifetime = pd.read_parquet(result.output_dir / "contact_lifetime.parquet")
    contact_map = np.load(result.output_dir / "contact_map.npy")
    summary = json.loads(result.summary_json.read_text(encoding="utf-8"))

    assert rg.columns.tolist() == ["frame", "time_ps", "rg"]
    assert ree.columns.tolist() == ["frame", "time_ps", "ree"]
    assert len(ps) == 2
    assert len(scaling) == 2
    assert msd["lag"].tolist() == [0, 1]
    assert lifetime["lag"].tolist() == [0, 1]
    assert contact_map.shape == (3, 3)
    assert summary["n_frames"] == 3
    assert summary["n_residues"] == 3
    assert summary["analysis"]["contact_cutoff_nm"] == 1.0
    assert "flory_fit" in summary


def test_analyze_trajectory_data_respects_optional_outputs(tmp_path) -> None:
    result = analyze_trajectory_data(
        _trajectory_data(),
        output_dir=tmp_path / "analysis",
        analysis_config=AnalysisConfig(observables=["rg", "ree", "contacts", "ps", "scaling"]),
    )

    assert "msd" not in result.outputs
    assert "contact_lifetime" not in result.outputs
    assert not (result.output_dir / "msd.parquet").exists()
    assert not (result.output_dir / "contact_lifetime.parquet").exists()
