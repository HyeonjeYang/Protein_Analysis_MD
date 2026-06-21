from __future__ import annotations

from idrptm.schema import SequenceConfig, WorkflowConfig
from idrptm.storage import estimate_dcd_bytes, estimate_project_storage


def test_100aa_1000_frames_estimates_about_1_2_mb() -> None:
    estimate = estimate_dcd_bytes(n_frames=1000, n_beads=100)

    assert estimate / 1_000_000 == 1.26


def test_project_storage_estimate_includes_runs() -> None:
    config = WorkflowConfig(
        project="storage_test",
        sequence=SequenceConfig(name="seq", sequence="A" * 10),
        calvados={"simulation": {"save_every_steps": 10, "n_frames": 5}},
    )

    estimate = estimate_project_storage(config)

    assert estimate.project == "storage_test"
    assert len(estimate.runs) == 1
    assert estimate.runs[0].n_beads == 10
    assert estimate.runs[0].n_frames == 5
