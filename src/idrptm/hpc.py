"""HPC scheduler scaffold placeholders."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SlurmConfig:
    """Minimal SLURM submission metadata."""

    job_name: str
    time: str = "24:00:00"
    partition: str | None = None
    gpus: int = 0
    cpus_per_task: int = 1
    memory: str = "4G"


def render_slurm_header(config: SlurmConfig) -> str:
    """Render a minimal SLURM header for Stage 2 expansion."""

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={config.job_name}",
        f"#SBATCH --time={config.time}",
        f"#SBATCH --cpus-per-task={config.cpus_per_task}",
        f"#SBATCH --mem={config.memory}",
    ]
    if config.partition:
        lines.append(f"#SBATCH --partition={config.partition}")
    if config.gpus:
        lines.append(f"#SBATCH --gres=gpu:{config.gpus}")
    return "\n".join(lines) + "\n"
