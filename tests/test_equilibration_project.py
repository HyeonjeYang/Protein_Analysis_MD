from __future__ import annotations

import json

import numpy as np
import pandas as pd

from idrptm.analysis.equilibration import equilibration_diagnostics, write_equilibration_outputs
from idrptm.project import clean_project, summarize_project_status


def test_equilibration_diagnostics_from_rg(tmp_path) -> None:
    diagnostics = equilibration_diagnostics(rg=np.array([1.0, 1.1, 1.0, 1.1]))

    assert diagnostics["available"] is True
    assert diagnostics["metric"] == "rg"
    assert "recommended_discard_frames" in diagnostics


def test_write_equilibration_outputs(tmp_path) -> None:
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir()
    pd.DataFrame({"frame": [0, 1, 2], "rg": [1.0, 1.1, 1.2]}).to_parquet(
        analysis_dir / "timeseries_rg.parquet",
        index=False,
    )

    outputs = write_equilibration_outputs(analysis_dir)

    assert outputs["equilibration_diagnostics"].exists()
    assert outputs["equilibration_blocks"].exists()


def test_project_status_and_clean(tmp_path) -> None:
    project = tmp_path / "project"
    run = project / "runs" / "seq__WT"
    analysis = run / "analysis"
    analysis.mkdir(parents=True)
    (run / "run.py").write_text("", encoding="utf-8")
    (run / "run_status.json").write_text(
        json.dumps({"status": "completed"}) + "\n",
        encoding="utf-8",
    )
    (analysis / "summary.json").write_text("{}\n", encoding="utf-8")
    (analysis / "cache_manifest.json").write_text("{}\n", encoding="utf-8")

    status = summarize_project_status(project)
    removed = clean_project(project, yes=True)

    assert status.counts["completed"] == 1
    assert status.counts["analyzed"] == 1
    assert analysis / "cache_manifest.json" in removed
