"""Run naming, parameter snapshots, and execution provenance helpers."""

from __future__ import annotations

import json
import os
import platform
import re
import socket
import sys
from datetime import datetime
from pathlib import Path
from pprint import pformat
from typing import Any


def slugify(value: object, *, default: str = "trajectory") -> str:
    """Return a filesystem-friendly identifier."""

    text = str(value or "").strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._-")
    return text or default


def timestamp_label(
    *,
    now: datetime | None = None,
    fmt: str = "%Y%m%d_%H%M%S",
) -> str:
    """Return a compact local timestamp label for run directories."""

    return (now or datetime.now()).strftime(fmt)


def parameter_snapshot(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Flatten nested parameters into a dict with value and type metadata."""

    flattened: dict[str, dict[str, Any]] = {}
    for key, value in _flatten(payload):
        flattened[key] = {
            "type": type(value).__name__,
            "value": _jsonable(value),
        }
    return flattened


def write_parameter_txt(
    path: str | Path,
    payload: dict[str, Any],
    *,
    title: str = "protein_analysis_md parameter snapshot",
) -> Path:
    """Write a paramdict-style text snapshot with explicit value types."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    paramdict = parameter_snapshot(payload)
    output.write_text(
        "\n".join(
            [
                f"# {title}",
                "# Format: paramdict[key] = {'type': python_type, 'value': value}",
                "paramdict = ",
                pformat(paramdict, sort_dicts=True, width=100),
                "",
            ]
        ),
        encoding="utf-8",
    )
    return output


def execution_environment() -> dict[str, Any]:
    """Collect lightweight local machine and Python execution metadata."""

    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "cwd": str(Path.cwd()),
        "pid": os.getpid(),
        "openmm_platforms": _openmm_platforms(),
    }


def build_trajectory_folder_name(
    *,
    project_name: str,
    protein_hint: str | None = None,
    preset: str | None = None,
    total_time_ns: float | None = None,
    replicates: int | None = None,
    ptm_mode: str | None = None,
    cleavage_mode: str | None = None,
    traj_name: str | None = None,
    traj_flag: str | None = None,
    include_timestamp: bool = True,
    timestamp_format: str = "%Y%m%d_%H%M%S",
    now: datetime | None = None,
) -> str:
    """Build a descriptive trajectory project directory name."""

    if traj_name:
        parts = [slugify(traj_name)]
    else:
        parts = [slugify(protein_hint or project_name)]
        if preset:
            parts.append(slugify(preset))
        if total_time_ns is not None:
            parts.append(slugify(f"{total_time_ns:g}ns"))
        if replicates and replicates > 1:
            parts.append(slugify(f"rep{replicates}"))
        if ptm_mode and ptm_mode not in {"none", "wt"}:
            parts.append(slugify(ptm_mode))
        if cleavage_mode and cleavage_mode != "none":
            parts.append(slugify(cleavage_mode))
    if include_timestamp:
        parts.append(timestamp_label(now=now, fmt=timestamp_format))
    if traj_flag:
        parts.append(slugify(traj_flag))
    return "__".join(part for part in parts if part)


def _flatten(
    payload: Any,
    *,
    prefix: str = "",
) -> list[tuple[str, Any]]:
    if isinstance(payload, dict):
        rows: list[tuple[str, Any]] = []
        for key in sorted(payload):
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(_flatten(payload[key], prefix=child_prefix))
        return rows
    if isinstance(payload, (list, tuple)):
        if all(not isinstance(item, dict | list | tuple) for item in payload):
            return [(prefix, payload)]
        rows = []
        for index, item in enumerate(payload):
            rows.extend(_flatten(item, prefix=f"{prefix}[{index}]"))
        return rows
    return [(prefix, payload)]


def _jsonable(value: Any) -> Any:
    try:
        json.dumps(value)
    except TypeError:
        if isinstance(value, Path):
            return str(value)
        return repr(value)
    return value


def _openmm_platforms() -> list[str]:
    try:
        from openmm import Platform
    except Exception:
        return []
    try:
        return [
            Platform.getPlatform(index).getName()
            for index in range(Platform.getNumPlatforms())
        ]
    except Exception:
        return []
