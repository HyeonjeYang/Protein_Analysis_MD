from __future__ import annotations

import csv
import json

from idrptm.visualization.pymol import export_pymol_project


def _write_project(tmp_path):
    project = tmp_path / "project"
    run_dir = project / "runs" / "toy__pSer2"
    run_dir.mkdir(parents=True)
    (run_dir / "top.pdb").write_text(
        "ATOM      1  CA  SER A   2       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (run_dir / "toy__pSer2.dcd").write_bytes(b"DCD")
    (run_dir / "metadata.json").write_text(
        json.dumps({"ptm_sites_1based": "2"}) + "\n",
        encoding="utf-8",
    )
    (run_dir / "parameters.txt").write_text("paramdict = {}\n", encoding="utf-8")
    with (project / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "metadata_path"])
        writer.writeheader()
        writer.writerow(
            {
                "variant_id": "toy__pSer2",
                "metadata_path": "runs/toy__pSer2/metadata.yaml",
            }
        )
    return project


def test_pymol_export_organizes_run_assets(tmp_path) -> None:
    project = _write_project(tmp_path)

    result = export_pymol_project(project, mode="copy")

    assert len(result.runs) == 1
    run = result.runs[0]
    assert run.status == "ready"
    assert (run.output_dir / "topology.pdb").read_text(encoding="utf-8").startswith("ATOM")
    assert (run.output_dir / "trajectory.dcd").read_bytes() == b"DCD"
    script = (run.output_dir / "load.pml").read_text(encoding="utf-8")
    assert "load topology.pdb" in script
    assert "load_traj trajectory.dcd" in script
    assert "select toy_pSer2_ptm_sites" in script
    assert result.manifest_csv.is_file()
    assert result.load_all_pml.is_file()
    assert result.readme.is_file()


def test_pymol_export_can_include_missing_runs(tmp_path) -> None:
    project = tmp_path / "project"
    run_dir = project / "runs" / "missing"
    run_dir.mkdir(parents=True)
    with (project / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "metadata_path"])
        writer.writeheader()
        writer.writerow({"variant_id": "missing", "metadata_path": "runs/missing/metadata.yaml"})

    result = export_pymol_project(project, include_missing=True)

    assert len(result.runs) == 1
    assert result.runs[0].status == "missing_topology_trajectory"
    assert (result.runs[0].output_dir / "load.pml").is_file()
