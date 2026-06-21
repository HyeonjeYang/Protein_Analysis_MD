"""Runtime environment diagnostics for local vs Remote-SSH execution."""

from __future__ import annotations

import getpass
import importlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
from pathlib import Path
from typing import Any

ENVIRONMENT_VARIABLES = (
    "SSH_CONNECTION",
    "SSH_CLIENT",
    "VSCODE_IPC_HOOK_CLI",
    "VSCODE_GIT_IPC_HANDLE",
    "TERM_PROGRAM",
    "CONDA_PREFIX",
    "VIRTUAL_ENV",
    "CUDA_VISIBLE_DEVICES",
)

REMOTE_INTERPRETATION = "This command appears to be running on a remote SSH host"
LOCAL_INTERPRETATION = "This command appears to be running locally"
AMBIGUOUS_INTERPRETATION = "Could not determine confidently"


def interpret_environment(env: dict[str, str] | None = None) -> str:
    """Return a conservative local-vs-remote interpretation."""

    values = dict(os.environ if env is None else env)
    if values.get("SSH_CONNECTION") or values.get("SSH_CLIENT"):
        return REMOTE_INTERPRETATION
    vscode_markers = values.get("VSCODE_IPC_HOOK_CLI") or values.get("VSCODE_GIT_IPC_HANDLE")
    if vscode_markers and values.get("TERM_PROGRAM") == "vscode":
        return AMBIGUOUS_INTERPRETATION
    return LOCAL_INTERPRETATION


def environment_appears_remote(env: dict[str, str] | None = None) -> bool:
    """Return true when SSH markers indicate remote execution."""

    return interpret_environment(env) == REMOTE_INTERPRETATION


def collect_environment_info(
    *,
    cwd: str | Path | None = None,
    env: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Collect diagnostic information without requiring optional tools."""

    root = Path(cwd or Path.cwd())
    env_values = dict(os.environ if env is None else env)
    selected_env = {name: env_values.get(name) for name in ENVIRONMENT_VARIABLES}
    openmm = _openmm_info()
    return {
        "hostname": socket.gethostname(),
        "fqdn": socket.getfqdn(),
        "cwd": str(root),
        "user": getpass.getuser(),
        "python_executable": sys.executable,
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "inside_vscode_remote_markers": _has_vscode_remote_markers(env_values),
        "environment_variables": selected_env,
        "cpu_count": os.cpu_count(),
        "memory": _memory_info(),
        "gpu": _gpu_info(),
        "calvados": _import_info("calvados"),
        "openmm": openmm,
        "openmm_platforms": openmm.get("platforms", []),
        "git": _git_info(root),
        "interpretation": interpret_environment(env_values),
    }


def format_environment_report(info: dict[str, Any]) -> str:
    """Format environment diagnostics for terminal output."""

    env = info.get("environment_variables", {})
    memory = info.get("memory") or {}
    gpu = info.get("gpu") or {}
    git = info.get("git") or {}
    lines = [
        f"hostname: {info.get('hostname')}",
        f"fqdn: {info.get('fqdn')}",
        f"cwd: {info.get('cwd')}",
        f"user: {info.get('user')}",
        f"python_executable: {info.get('python_executable')}",
        f"python_version: {info.get('python_version')}",
        f"platform: {info.get('platform')}",
        f"system: {info.get('system')} {info.get('release')} {info.get('machine')}",
        f"inside_vscode_remote_markers: {info.get('inside_vscode_remote_markers')}",
        "environment_variables:",
    ]
    lines.extend(f"  {name}: {env.get(name)}" for name in ENVIRONMENT_VARIABLES)
    lines.extend(
        [
            f"cpu_count: {info.get('cpu_count')}",
            f"memory_total_bytes: {memory.get('total')}",
            f"memory_available_bytes: {memory.get('available')}",
            f"gpu_available: {gpu.get('available')}",
            f"gpu_summary: {gpu.get('summary')}",
            f"calvados_importable: {info.get('calvados', {}).get('importable')}",
            f"calvados_error: {info.get('calvados', {}).get('error')}",
            f"openmm_importable: {info.get('openmm', {}).get('importable')}",
            f"openmm_error: {info.get('openmm', {}).get('error')}",
            f"openmm_platforms: {', '.join(info.get('openmm_platforms') or []) or None}",
            f"git_branch: {git.get('branch')}",
            f"git_commit: {git.get('commit')}",
            f"interpretation: {info.get('interpretation')}",
        ]
    )
    return "\n".join(lines)


def write_environment_json(info: dict[str, Any], path: str | Path) -> Path:
    """Write environment diagnostics as JSON."""

    output_path = Path(path)
    output_path.write_text(json.dumps(info, indent=2, default=str) + "\n", encoding="utf-8")
    return output_path


def _has_vscode_remote_markers(env: dict[str, str]) -> bool:
    return bool(env.get("SSH_CONNECTION") or env.get("SSH_CLIENT")) and bool(
        env.get("VSCODE_IPC_HOOK_CLI") or env.get("VSCODE_GIT_IPC_HANDLE")
    )


def _memory_info() -> dict[str, int] | None:
    try:
        psutil = importlib.import_module("psutil")
    except Exception:
        return None
    try:
        memory = psutil.virtual_memory()
    except Exception:
        return None
    return {"total": int(memory.total), "available": int(memory.available)}


def _gpu_info() -> dict[str, Any]:
    executable = shutil.which("nvidia-smi")
    if executable is None:
        return {"available": False, "summary": "nvidia-smi not found"}
    try:
        completed = subprocess.run(
            [
                executable,
                "--query-gpu=name,memory.total,memory.free",
                "--format=csv,noheader",
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception as exc:
        return {"available": False, "summary": f"nvidia-smi failed: {exc}"}
    if completed.returncode != 0:
        return {"available": False, "summary": completed.stderr.strip() or "nvidia-smi failed"}
    return {"available": True, "summary": completed.stdout.strip()}


def _import_info(module_name: str) -> dict[str, Any]:
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        return {"importable": False, "module": module_name, "error": str(exc)}
    return {
        "importable": True,
        "module": module_name,
        "version": getattr(module, "__version__", None),
        "error": None,
    }


def _openmm_info() -> dict[str, Any]:
    info = _import_info("openmm")
    if not info["importable"]:
        return info | {"platforms": []}
    try:
        openmm = importlib.import_module("openmm")
        platforms = [
            openmm.Platform.getPlatform(index).getName()
            for index in range(openmm.Platform.getNumPlatforms())
        ]
    except Exception as exc:
        return info | {"platforms": [], "error": str(exc)}
    return info | {"platforms": platforms}


def _git_info(cwd: Path) -> dict[str, str | None]:
    return {
        "branch": _git_output(cwd, "rev-parse", "--abbrev-ref", "HEAD"),
        "commit": _git_output(cwd, "rev-parse", "--short", "HEAD"),
    }


def _git_output(cwd: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=cwd,
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None
