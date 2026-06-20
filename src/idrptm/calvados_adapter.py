"""Adapter boundaries for external CALVADOS runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CalvadosRunDirectory:
    """Paths that identify a prepared CALVADOS run directory."""

    path: Path
    config_yaml: Path
    components_yaml: Path
    run_script: Path


def prepare_calvados_run_directory(output_dir: str | Path) -> CalvadosRunDirectory:
    """Declare where Stage 2 will generate CALVADOS inputs.

    This function intentionally does not create full CALVADOS inputs yet.
    """

    run_dir = Path(output_dir)
    return CalvadosRunDirectory(
        path=run_dir,
        config_yaml=run_dir / "config.yaml",
        components_yaml=run_dir / "components.yaml",
        run_script=run_dir / "run.py",
    )
