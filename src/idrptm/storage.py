"""Trajectory storage estimation utilities."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from idrptm.design import DesignedVariant, design_variants
from idrptm.schema import WorkflowConfig, load_config

DCD_OVERHEAD_FACTOR = 1.05
BYTES_PER_FLOAT32 = 4
COORDINATES_PER_BEAD = 3


@dataclass(frozen=True)
class RunStorageEstimate:
    """Storage estimate for one designed run."""

    run_id: str
    n_beads: int
    n_frames: int
    dcd_bytes: float

    @property
    def dcd_mb(self) -> float:
        return self.dcd_bytes / 1_000_000.0


@dataclass(frozen=True)
class ProjectStorageEstimate:
    """Storage estimate for all designed runs in a project."""

    project: str
    runs: tuple[RunStorageEstimate, ...]
    warn_gb: float = 10.0
    strong_warn_gb: float = 100.0

    @property
    def total_dcd_bytes(self) -> float:
        return sum(run.dcd_bytes for run in self.runs)

    @property
    def total_dcd_gb(self) -> float:
        return self.total_dcd_bytes / 1_000_000_000.0

    @property
    def warning_level(self) -> str:
        if self.total_dcd_gb >= self.strong_warn_gb:
            return "strong"
        if self.total_dcd_gb >= self.warn_gb:
            return "warn"
        return "ok"

    def to_dict(self) -> dict[str, object]:
        return {
            "project": self.project,
            "canonical_units": {
                "storage": "bytes",
                "storage_report": "MB/GB",
            },
            "total_dcd_bytes": self.total_dcd_bytes,
            "total_dcd_gb": self.total_dcd_gb,
            "warning_level": self.warning_level,
            "runs": [
                {
                    "run_id": run.run_id,
                    "n_beads": run.n_beads,
                    "n_frames": run.n_frames,
                    "dcd_bytes": run.dcd_bytes,
                    "dcd_mb": run.dcd_mb,
                }
                for run in self.runs
            ],
        }


def estimate_dcd_bytes(n_frames: int, n_beads: int) -> float:
    """Estimate DCD coordinate storage bytes for float32 XYZ coordinates."""

    return n_frames * n_beads * COORDINATES_PER_BEAD * BYTES_PER_FLOAT32 * DCD_OVERHEAD_FACTOR


def estimate_project_storage(config: WorkflowConfig) -> ProjectStorageEstimate:
    """Estimate trajectory storage for each designed run in a workflow."""

    n_frames = int(config.calvados.simulation.n_frames or 0)
    runs = tuple(
        RunStorageEstimate(
            run_id=variant.variant_id,
            n_beads=_variant_bead_count(variant),
            n_frames=n_frames,
            dcd_bytes=estimate_dcd_bytes(n_frames, _variant_bead_count(variant)),
        )
        for variant in design_variants(config)
    )
    return ProjectStorageEstimate(project=config.project, runs=runs)


def estimate_from_config_file(
    config_path: str | Path,
    *,
    output_dir: str | Path | None = None,
) -> ProjectStorageEstimate:
    """Load a config, estimate storage, and write ``storage_estimate.json``."""

    config = load_config(config_path)
    estimate = estimate_project_storage(config)
    root = Path(output_dir) if output_dir is not None else config.runner.work_dir
    root.mkdir(parents=True, exist_ok=True)
    write_storage_estimate(estimate, root / "storage_estimate.json")
    return estimate


def write_storage_estimate(estimate: ProjectStorageEstimate, path: str | Path) -> Path:
    """Write a storage estimate JSON file."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(estimate.to_dict(), indent=2) + "\n", encoding="utf-8")
    return output_path


def format_storage_table(estimate: ProjectStorageEstimate) -> str:
    """Return a readable storage-estimate table."""

    lines = [
        "run_id,n_beads,n_frames,dcd_mb",
        *[
            f"{run.run_id},{run.n_beads},{run.n_frames},{run.dcd_mb:.3f}"
            for run in estimate.runs
        ],
        f"TOTAL,,,{estimate.total_dcd_gb * 1000:.3f}",
        f"warning_level,{estimate.warning_level},,",
    ]
    return "\n".join(lines)


def _variant_bead_count(variant: DesignedVariant) -> int:
    return sum(
        len(component.simulation_sequence) * component.copies
        for component in variant.components
    )
