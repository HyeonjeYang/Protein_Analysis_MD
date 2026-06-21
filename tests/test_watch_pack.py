from __future__ import annotations

import csv
import json
import tarfile

from idrptm.bundle import pack_project
from idrptm.watch import format_watch, summarize_watch


def _project(tmp_path):
    project = tmp_path / "project"
    run = project / "runs" / "toy__WT"
    run.mkdir(parents=True)
    with (project / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "metadata_path"])
        writer.writeheader()
        writer.writerow({"variant_id": "toy__WT", "metadata_path": "runs/toy__WT/metadata.yaml"})
    (project / "project.lock.yaml").write_text("project: toy\n", encoding="utf-8")
    (run / "run.py").write_text("print('run')\n", encoding="utf-8")
    (run / "run_status.json").write_text(
        json.dumps(
            {
                "status": "submitted",
                "progress": {"frames_written_estimate": 2},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (run / "config.yaml").write_text("steps: 100\nwfreq: 10\n", encoding="utf-8")
    (run / "metadata.json").write_text("{}\n", encoding="utf-8")
    (run / "parameters.txt").write_text("paramdict = {}\n", encoding="utf-8")
    (run / "top.pdb").write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        encoding="utf-8",
    )
    (run / "toy__WT.dcd").write_bytes(b"DCD")
    return project


def test_watch_formats_project_progress(tmp_path) -> None:
    project = _project(tmp_path)

    snapshot = summarize_watch(project)
    text = format_watch(snapshot)

    assert snapshot.n_runs == 1
    assert snapshot.ready_trajectories == 1
    assert snapshot.frame_progress == (2, 10)
    assert "frames_estimate: 2/10" in text


def test_pack_excludes_trajectories_by_default(tmp_path) -> None:
    project = _project(tmp_path)

    result = pack_project(project)

    with tarfile.open(result.archive_path, "r:gz") as tar:
        names = tar.getnames()
    assert any(name.endswith("parameters.txt") for name in names)
    assert not any(name.endswith(".dcd") for name in names)
    assert not any(name.endswith("top.pdb") for name in names)


def test_pack_can_include_trajectories(tmp_path) -> None:
    project = _project(tmp_path)

    result = pack_project(project, output=tmp_path / "with_traj.tar.gz", include_trajectories=True)

    with tarfile.open(result.archive_path, "r:gz") as tar:
        names = tar.getnames()
    assert any(name.endswith(".dcd") for name in names)
    assert any(name.endswith("top.pdb") for name in names)
