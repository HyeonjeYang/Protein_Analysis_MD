"""Local detached execution helpers with hardware-aware backend selection."""

from __future__ import annotations

import csv
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import yaml

BackendChoice = Literal["auto", "CPU", "CUDA", "OpenCL"]
ResolvedBackend = Literal["CPU", "CUDA", "OpenCL"]
TerminalChoice = Literal["auto", "tmux", "byobu", "none"]
ResolvedTerminal = Literal["tmux", "byobu", "none"]


@dataclass(frozen=True)
class HardwareInfo:
    """Small hardware summary for local launch decisions."""

    system: str
    machine: str
    cpu_count: int
    performance_cores: int | None
    chip: str | None
    gpu_summary: str | None
    openmm_platforms: tuple[str, ...]


@dataclass(frozen=True)
class LaunchRecommendation:
    """Recommended local execution settings."""

    backend: ResolvedBackend
    simulation_parallel: int
    analysis_parallel: int
    terminal: ResolvedTerminal
    cpu_threads_per_run: int
    reason: str
    hardware: HardwareInfo

    def summary(self) -> str:
        """Return a compact human-readable recommendation."""

        platforms = ", ".join(self.hardware.openmm_platforms) or "unknown"
        return "\n".join(
            [
                "Recommended local execution:",
                f"  backend: {self.backend}",
                f"  simulation_parallel: {self.simulation_parallel}",
                f"  analysis_parallel: {self.analysis_parallel}",
                f"  terminal: {self.terminal}",
                f"  cpu_threads_per_run: {self.cpu_threads_per_run}",
                f"  OpenMM platforms: {platforms}",
                f"  hardware: {self.hardware.chip or self.hardware.machine}",
                f"  reason: {self.reason}",
            ]
        )


def collect_hardware_info() -> HardwareInfo:
    """Collect local hardware and OpenMM platform details without failing hard."""

    system = platform.system()
    machine = platform.machine()
    cpu_count = os.cpu_count() or 1
    chip = None
    gpu_summary = None
    performance_cores = None
    if system == "Darwin":
        profiler = _run_text(["system_profiler", "SPHardwareDataType", "SPDisplaysDataType"])
        chip_match = re.search(r"Chip:\s*(.+)", profiler)
        if chip_match:
            chip = chip_match.group(1).strip()
        core_match = re.search(
            r"Total Number of Cores:\s*\d+\s*\((\d+)\s+Performance",
            profiler,
        )
        if core_match:
            performance_cores = int(core_match.group(1))
        gpu_match = re.search(r"Chipset Model:\s*(.+)", profiler)
        if gpu_match:
            gpu_summary = gpu_match.group(1).strip()
    return HardwareInfo(
        system=system,
        machine=machine,
        cpu_count=cpu_count,
        performance_cores=performance_cores,
        chip=chip,
        gpu_summary=gpu_summary,
        openmm_platforms=_openmm_platforms(),
    )


def recommend_local_execution(
    *,
    backend: BackendChoice = "auto",
    simulation_parallel: str | int = "auto",
    analysis_parallel: str | int = "auto",
    terminal: TerminalChoice = "auto",
    hardware: HardwareInfo | None = None,
) -> LaunchRecommendation:
    """Choose backend and parallelism from requested settings and hardware."""

    info = hardware or collect_hardware_info()
    selected_backend, reason = _select_backend(backend, info)
    sim_parallel = _resolve_parallel(
        simulation_parallel,
        backend=selected_backend,
        hardware=info,
        for_analysis=False,
    )
    ana_parallel = _resolve_parallel(
        analysis_parallel,
        backend=selected_backend,
        hardware=info,
        for_analysis=True,
    )
    return LaunchRecommendation(
        backend=selected_backend,
        simulation_parallel=sim_parallel,
        analysis_parallel=ana_parallel,
        terminal=_select_terminal(terminal),
        cpu_threads_per_run=1,
        reason=reason,
        hardware=info,
    )


def apply_backend_to_project(
    project_dir: str | Path,
    *,
    backend: ResolvedBackend,
    cpu_threads_per_run: int = 1,
) -> tuple[Path, ...]:
    """Update prepared run config.yaml files to use the selected backend."""

    updated: list[Path] = []
    for run_dir in _run_dirs(Path(project_dir)):
        config_path = run_dir / "config.yaml"
        if not config_path.exists():
            continue
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        payload["platform"] = backend
        if backend == "CPU":
            payload["threads"] = cpu_threads_per_run
        config_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
        updated.append(config_path)
    return tuple(updated)


def launch_local_pipeline(
    project_dir: str | Path,
    *,
    recommendation: LaunchRecommendation,
    python_executable: str,
    session_name: str | None = None,
    replace_existing: bool = True,
) -> dict[str, object]:
    """Launch a detached local pipeline in tmux/byobu or a plain session."""

    root = Path(project_dir).resolve()
    repo_root = _repo_root_for_project(root)
    if replace_existing:
        stop_existing_pipeline(root, session_name=session_name)
    session = sanitize_session_name(session_name or f"pamd_{root.name}")
    command = _worker_command(
        project_dir=root,
        python_executable=python_executable,
        simulation_parallel=recommendation.simulation_parallel,
        analysis_parallel=recommendation.analysis_parallel,
    )
    env = _worker_env(repo_root)
    if recommendation.terminal == "tmux":
        subprocess.run(
            ["tmux", "new-session", "-d", "-s", session, *command],
            check=True,
            cwd=repo_root,
            env=env,
        )
        pid = None
    elif recommendation.terminal == "byobu":
        subprocess.run(
            ["byobu", "new-session", "-d", "-s", session, *command],
            check=True,
            cwd=repo_root,
            env=env,
        )
        pid = None
    else:
        log_path = root / "pipeline.nohup.log"
        log = log_path.open("ab", buffering=0)
        process = subprocess.Popen(
            command,
            cwd=repo_root,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            start_new_session=True,
        )
        (root / "pipeline.pid").write_text(f"{process.pid}\n", encoding="utf-8")
        pid = process.pid
    (root / "launch_settings.json").write_text(
        json.dumps(
            {
                "backend": recommendation.backend,
                "simulation_parallel": recommendation.simulation_parallel,
                "analysis_parallel": recommendation.analysis_parallel,
                "terminal": recommendation.terminal,
                "session_name": session,
                "pid": pid,
                "launched_at": _now(),
                "repo_root": str(repo_root),
                "hardware": recommendation.hardware.__dict__,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"session_name": session, "pid": pid, "command": command}


def run_pipeline_worker(
    project_dir: str | Path,
    *,
    simulation_parallel: int,
    analysis_parallel: int,
    python_executable: str,
) -> int:
    """Run simulations, analyses, comparison, and report generation."""

    root = Path(project_dir).resolve()
    repo_root = _repo_root_for_project(root)
    env = _worker_env(repo_root)
    logs = root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    status_path = root / "pipeline_status.json"
    run_dirs = _run_dirs(root)
    _write_status(
        status_path,
        {
            "stage": "simulation",
            "status": "running",
            "started_at": _now(),
            "project_dir": str(root),
            "n_runs": len(run_dirs),
            "simulation_parallel": simulation_parallel,
            "analysis_parallel": analysis_parallel,
        },
    )
    simulation_results = _run_pool(
        [
            (
                run_dir.name,
                [
                    *_cli_command(python_executable),
                    "run",
                    str(run_dir),
                    "--phase",
                    "all",
                    "--execute",
                    "--python",
                    python_executable,
                ],
                logs / f"{run_dir.name}.simulation.log",
            )
            for run_dir in run_dirs
        ],
        max_parallel=simulation_parallel,
        status_path=status_path,
        stage="simulation",
        cwd=repo_root,
        env=env,
    )
    failed = [item for item in simulation_results if item["returncode"] != 0]
    if failed:
        _write_status(status_path, _failure_payload("simulation", failed, simulation_results))
        return 1

    lock = root / "project.lock.yaml"
    analysis_results = _run_pool(
        [
            (
                run_dir.name,
                [
                    *_cli_command(python_executable),
                    "analyze",
                    str(run_dir),
                    "--config",
                    str(lock),
                    "--force",
                ],
                logs / f"{run_dir.name}.analysis.log",
            )
            for run_dir in run_dirs
        ],
        max_parallel=analysis_parallel,
        status_path=status_path,
        stage="analysis",
        cwd=repo_root,
        env=env,
    )
    failed = [item for item in analysis_results if item["returncode"] != 0]
    if failed:
        _write_status(status_path, _failure_payload("analysis", failed, analysis_results))
        return 1

    final_results = _run_pool(
        _final_commands(
            root=root,
            logs=logs,
            python_executable=python_executable,
            visualization=_visualization_enabled(root),
        ),
        max_parallel=1,
        status_path=status_path,
        stage="report",
        cwd=repo_root,
        env=env,
    )
    failed = [item for item in final_results if item["returncode"] != 0]
    _write_status(
        status_path,
        {
            "stage": "done" if not failed else "report",
            "status": "completed" if not failed else "failed",
            "ended_at": _now(),
            "simulation": simulation_results,
            "analysis": analysis_results,
            "final": final_results,
        },
    )
    return 1 if failed else 0


def stop_existing_pipeline(project_dir: str | Path, *, session_name: str | None = None) -> None:
    """Best-effort stop for a previously launched local pipeline."""

    root = Path(project_dir).resolve()
    session = sanitize_session_name(session_name or f"pamd_{root.name}")
    if shutil.which("tmux"):
        subprocess.run(["tmux", "kill-session", "-t", session], check=False)
    if shutil.which("byobu"):
        subprocess.run(["byobu", "kill-session", "-t", session], check=False)
    pid_path = root / "pipeline.pid"
    if pid_path.exists():
        try:
            pid = int(pid_path.read_text(encoding="utf-8").strip())
            os.killpg(pid, signal.SIGTERM)
        except Exception:
            try:
                os.kill(pid, signal.SIGTERM)
            except Exception:
                pass
        pid_path.unlink(missing_ok=True)


def clean_interrupted_outputs(project_dir: str | Path, run_names: tuple[str, ...] = ()) -> None:
    """Remove partial trajectory/checkpoint outputs from interrupted runs."""

    patterns = (
        "*.dcd",
        "*.xml",
        "top.pdb",
        "bonds_*.txt",
        "*.log",
        "run_status.json",
        "restart.chk",
    )
    run_dirs = _run_dirs(Path(project_dir))
    selected = {
        name for name in run_names
    }
    for run_dir in run_dirs:
        if selected and run_dir.name not in selected:
            continue
        for pattern in patterns:
            for path in run_dir.glob(pattern):
                path.unlink(missing_ok=True)


def sanitize_session_name(name: str) -> str:
    """Return a tmux/byobu-safe session name."""

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "pamd"


def _select_backend(
    requested: BackendChoice,
    hardware: HardwareInfo,
) -> tuple[ResolvedBackend, str]:
    platforms = set(hardware.openmm_platforms)
    if requested != "auto":
        if requested not in platforms:
            raise ValueError(
                f"Requested OpenMM backend {requested!r} is not available: "
                f"{sorted(platforms)}"
            )
        return requested, "Backend explicitly requested by user."
    if "CUDA" in platforms:
        return "CUDA", "CUDA is available and is normally best for production MD."
    if hardware.system == "Darwin":
        return (
            "CPU",
            "macOS/Apple Silicon has no CUDA; OpenCL may be listed but is often unavailable "
            "for CALVADOS/OpenMM contexts, so CPU is the conservative automatic choice.",
        )
    if "OpenCL" in platforms and "CPU" not in platforms:
        return "OpenCL", "OpenCL is the only accelerator-style backend detected."
    return "CPU", "CPU is available and no reliable CUDA backend was detected."


def _resolve_parallel(
    value: str | int,
    *,
    backend: ResolvedBackend,
    hardware: HardwareInfo,
    for_analysis: bool,
) -> int:
    if value != "auto":
        return max(1, int(value))
    if for_analysis:
        return max(1, min(4, hardware.performance_cores or hardware.cpu_count))
    if backend == "CPU":
        cores = hardware.performance_cores or max(1, hardware.cpu_count - 2)
        return max(1, cores)
    return 1


def _select_terminal(requested: TerminalChoice) -> ResolvedTerminal:
    if requested == "auto":
        if shutil.which("byobu"):
            return "byobu"
        if shutil.which("tmux"):
            return "tmux"
        return "none"
    if requested != "none" and shutil.which(requested) is None:
        raise ValueError(f"Requested terminal {requested!r} is not installed.")
    return requested


def _openmm_platforms() -> tuple[str, ...]:
    try:
        from openmm import Platform
    except Exception:
        return ()
    return tuple(
        Platform.getPlatform(index).getName() for index in range(Platform.getNumPlatforms())
    )


def _run_text(command: list[str]) -> str:
    try:
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
    except Exception:
        return ""
    return completed.stdout


def _run_dirs(root: Path) -> tuple[Path, ...]:
    manifest = root / "manifest.csv"
    if manifest.exists():
        dirs: list[Path] = []
        with manifest.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("metadata_path"):
                    dirs.append((root / row["metadata_path"]).parent)
        return tuple(sorted(dict.fromkeys(dirs)))
    runs = root / "runs"
    if runs.exists():
        return tuple(sorted(path for path in runs.iterdir() if (path / "run.py").exists()))
    return ()


def _worker_command(
    *,
    project_dir: Path,
    python_executable: str,
    simulation_parallel: int,
    analysis_parallel: int,
) -> list[str]:
    return [
        *_cli_command(python_executable),
        "local-pipeline-worker",
        str(project_dir),
        "--simulation-parallel",
        str(simulation_parallel),
        "--analysis-parallel",
        str(analysis_parallel),
        "--python",
        python_executable,
    ]


def _final_commands(
    *,
    root: Path,
    logs: Path,
    python_executable: str,
    visualization: bool,
) -> list[tuple[str, list[str], Path]]:
    commands = [
        (
            "compare",
            [*_cli_command(python_executable), "compare", str(root)],
            logs / "compare.log",
        )
    ]
    if visualization:
        commands.extend(
            [
                (
                    "report",
                    [*_cli_command(python_executable), "report", str(root)],
                    logs / "report.log",
                ),
                (
                    "pymol",
                    [
                        *_cli_command(python_executable),
                        "pymol",
                        str(root),
                        "--force",
                    ],
                    logs / "pymol.log",
                ),
            ]
        )
    commands.append(
        (
            "dashboard",
            [*_cli_command(python_executable), "dashboard", str(root)],
            logs / "dashboard.log",
        )
    )
    return commands


def _visualization_enabled(root: Path) -> bool:
    lock = root / "project.lock.yaml"
    if not lock.exists():
        return True
    try:
        payload = yaml.safe_load(lock.read_text(encoding="utf-8")) or {}
    except Exception:
        return True
    return bool(payload.get("visualization", True))


def _cli_command(python_executable: str) -> list[str]:
    return [python_executable, "-m", "idrptm.cli"]


def _repo_root_for_project(project_dir: Path) -> Path:
    """Find the repository root used for editable imports and relative files."""

    path = project_dir.resolve()
    for candidate in [path, *path.parents, Path.cwd().resolve(), *Path.cwd().resolve().parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "src").is_dir():
            return candidate
    return Path.cwd().resolve()


def _worker_env(repo_root: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(repo_root / ".mplconfig")
    env["PYTHONUNBUFFERED"] = "1"
    env["OPENMM_CPU_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def _run_pool(
    commands: list[tuple[str, list[str], Path]],
    *,
    max_parallel: int,
    status_path: Path,
    stage: str,
    cwd: Path,
    env: dict[str, str],
) -> list[dict[str, object]]:
    pending = list(commands)
    running: list[tuple[str, subprocess.Popen[bytes], object, Path, float]] = []
    completed: list[dict[str, object]] = []
    while pending or running:
        while pending and len(running) < max_parallel:
            name, command, log_path = pending.pop(0)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handle = log_path.open("ab")
            handle.write(f"\n[{_now()}] START {' '.join(command)}\n".encode())
            handle.flush()
            running.append(
                (
                    name,
                    subprocess.Popen(
                        command,
                        stdout=handle,
                        stderr=subprocess.STDOUT,
                        cwd=cwd,
                        env=env,
                    ),
                    handle,
                    log_path,
                    time.time(),
                )
            )
        next_running: list[tuple[str, subprocess.Popen[bytes], object, Path, float]] = []
        for name, process, handle, log_path, started_at in running:
            returncode = process.poll()
            if returncode is None:
                next_running.append((name, process, handle, log_path, started_at))
                continue
            elapsed_s = round(time.time() - started_at, 1)
            handle.write(f"[{_now()}] END returncode={returncode} elapsed_s={elapsed_s}\n".encode())
            handle.close()
            completed.append(
                {
                    "name": name,
                    "returncode": returncode,
                    "elapsed_s": elapsed_s,
                    "log": str(log_path),
                }
            )
        running = next_running
        _write_status(
            status_path,
            {
                "stage": stage,
                "status": "running",
                "updated_at": _now(),
                "pending": len(pending),
                "running": [name for name, *_ in running],
                "completed": completed,
            },
        )
        if pending or running:
            time.sleep(5)
    return completed


def _failure_payload(
    stage: str,
    failed: list[dict[str, object]],
    completed: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "stage": stage,
        "status": "failed",
        "ended_at": _now(),
        "failed": failed,
        "completed": completed,
    }


def _write_status(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()
