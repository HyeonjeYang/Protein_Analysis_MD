from __future__ import annotations

import csv

import yaml

from idrptm.local_pipeline import (
    HardwareInfo,
    _final_commands,
    apply_backend_to_project,
    recommend_local_execution,
)


def _hardware(
    *,
    system: str = "Darwin",
    platforms: tuple[str, ...] = ("Reference", "CPU", "OpenCL"),
    performance_cores: int | None = 6,
) -> HardwareInfo:
    return HardwareInfo(
        system=system,
        machine="arm64",
        cpu_count=8,
        performance_cores=performance_cores,
        chip="Apple M1 Pro",
        gpu_summary="Apple M1 Pro",
        openmm_platforms=platforms,
    )


def test_auto_backend_prefers_cpu_on_apple_silicon() -> None:
    recommendation = recommend_local_execution(
        hardware=_hardware(),
        backend="auto",
        terminal="none",
    )

    assert recommendation.backend == "CPU"
    assert recommendation.simulation_parallel == 6
    assert "Apple Silicon" in recommendation.reason


def test_auto_backend_prefers_cuda_when_available() -> None:
    recommendation = recommend_local_execution(
        hardware=_hardware(system="Linux", platforms=("Reference", "CPU", "CUDA")),
        backend="auto",
        terminal="none",
    )

    assert recommendation.backend == "CUDA"
    assert recommendation.simulation_parallel == 1


def test_explicit_opencl_backend_can_be_requested() -> None:
    recommendation = recommend_local_execution(
        hardware=_hardware(system="Linux", platforms=("Reference", "CPU", "OpenCL")),
        backend="OpenCL",
        terminal="none",
    )

    assert recommendation.backend == "OpenCL"


def test_apply_backend_updates_prepared_run_configs(tmp_path) -> None:
    project = tmp_path / "project"
    run_dir = project / "runs" / "rep001"
    run_dir.mkdir(parents=True)
    (run_dir / "metadata.json").write_text("{}\n", encoding="utf-8")
    (run_dir / "config.yaml").write_text(
        yaml.safe_dump({"platform": "CPU", "threads": 4}),
        encoding="utf-8",
    )
    project.mkdir(exist_ok=True)
    with (project / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id", "metadata_path"])
        writer.writeheader()
        writer.writerow({"run_id": "rep001", "metadata_path": "runs/rep001/metadata.json"})

    updated = apply_backend_to_project(project, backend="CPU", cpu_threads_per_run=1)
    payload = yaml.safe_load((run_dir / "config.yaml").read_text(encoding="utf-8"))

    assert updated == (run_dir / "config.yaml",)
    assert payload["platform"] == "CPU"
    assert payload["threads"] == 1


def test_final_commands_respect_visualization_flag(tmp_path) -> None:
    logs = tmp_path / "logs"

    without_visuals = _final_commands(
        root=tmp_path,
        logs=logs,
        python_executable="python",
        visualization=False,
    )
    with_visuals = _final_commands(
        root=tmp_path,
        logs=logs,
        python_executable="python",
        visualization=True,
    )

    assert [item[0] for item in without_visuals] == ["compare", "dashboard"]
    assert [item[0] for item in with_visuals] == ["compare", "report", "pymol", "dashboard"]
