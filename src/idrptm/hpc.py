"""HPC scheduler script generation helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from idrptm.runner import RunPhase, discover_run_directories


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
    """Render a minimal SLURM header."""

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


def write_slurm_array_script(
    project_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    phase: RunPhase = "all",
    python_executable: str = "python",
    config: SlurmConfig | None = None,
) -> Path:
    """Write a SLURM array script for all prepared runs in a project."""

    root = Path(project_dir)
    destination = Path(output_dir) if output_dir is not None else root
    destination.mkdir(parents=True, exist_ok=True)
    run_dirs = discover_run_directories(root, all_runs=True)
    if not run_dirs:
        raise ValueError(f"No prepared run directories found under {root}.")

    run_list = destination / "run_dirs.txt"
    run_list.write_text(
        "\n".join(str(path.resolve()) for path in run_dirs) + "\n",
        encoding="utf-8",
    )

    slurm_config = config or SlurmConfig(job_name=root.name or "protein_analysis_md")
    header = render_slurm_header(slurm_config).rstrip()
    script = destination / "run_slurm_array.sh"
    script.write_text(
        "\n".join(
            [
                header,
                f"#SBATCH --array=0-{len(run_dirs) - 1}",
                "",
                "set -euo pipefail",
                f"RUN_LIST={run_list.resolve()}",
                'RUN_DIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RUN_LIST")',
                f'pamd run "$RUN_DIR" --phase {phase} --execute --python "{python_executable}"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    return script
