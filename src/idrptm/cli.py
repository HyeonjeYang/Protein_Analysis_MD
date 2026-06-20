"""Command-line interface for idr-ptm-md."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from idrptm import __version__
from idrptm.calvados_adapter import prepare_from_config_file
from idrptm.design import design_from_config_file

app = typer.Typer(
    name="idrptm",
    help="Prepare, run, and analyze CALVADOS-backed IDR/PTM MD workflows.",
    no_args_is_help=True,
)


ConfigOption = Annotated[
    Path | None,
    typer.Option("--config", "-c", help="Workflow YAML configuration file."),
]


@app.callback()
def callback(
    version: Annotated[
        bool,
        typer.Option("--version", help="Show the idr-ptm-md version and exit."),
    ] = False,
) -> None:
    """Top-level CLI callback."""

    if version:
        typer.echo(f"idr-ptm-md {__version__}")
        raise typer.Exit()


@app.command("init")
def init_command(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Directory where a starter project would be created."),
    ] = Path("."),
    force: Annotated[
        bool,
        typer.Option("--force", help="Allow overwriting existing starter files."),
    ] = False,
) -> None:
    """Placeholder for creating a starter workflow directory."""

    typer.echo(
        f"Stage 1 placeholder: would initialize an idr-ptm-md project at {output} "
        f"(force={force})."
    )


@app.command("design")
def design_command(
    config: Annotated[
        Path,
        typer.Argument(help="Workflow YAML configuration file."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for design outputs."),
    ] = None,
) -> None:
    """Generate WT/PTM design manifest, FASTA files, and metadata stubs."""

    try:
        result = design_from_config_file(config, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Design failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote manifest: {result.manifest_path}")
    typer.echo(f"Wrote {len(result.fasta_paths)} variant FASTA file(s).")
    typer.echo(f"Wrote {len(result.metadata_paths)} per-run metadata stub(s).")


@app.command("prepare")
def prepare_command(
    config: Annotated[
        Path,
        typer.Argument(help="Workflow YAML configuration file."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="CALVADOS run-directory root."),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Print planned run directories without writing files."),
    ] = False,
) -> None:
    """Prepare CALVADOS run directories from a workflow config."""

    try:
        result = prepare_from_config_file(config, output_dir=output_dir, dry_run=dry_run)
    except Exception as exc:
        typer.echo(f"Prepare failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    if result.dry_run:
        typer.echo(f"Dry run: would prepare {len(result.run_directories)} run directories.")
        for run_dir in result.run_directories:
            typer.echo(f"Would prepare: {run_dir.path}")
        return

    typer.echo(f"Prepared {len(result.run_directories)} CALVADOS run directories.")
    typer.echo(f"Manifest: {result.manifest_path}")


@app.command("run")
def run_command(
    config: ConfigOption = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run/--execute",
            help="Print planned commands instead of executing them.",
        ),
    ] = True,
) -> None:
    """Placeholder for local/HPC execution."""

    typer.echo(f"Stage 1 placeholder: would run workflow from {config} (dry_run={dry_run}).")


@app.command("analyze")
def analyze_command(
    config: ConfigOption = None,
    trajectories: Annotated[
        list[Path] | None,
        typer.Option("--trajectory", "-t", help="Trajectory file to analyze."),
    ] = None,
) -> None:
    """Placeholder for trajectory analysis."""

    count = 0 if trajectories is None else len(trajectories)
    typer.echo(
        f"Stage 1 placeholder: would analyze {count} trajectory file(s) using config {config}."
    )


@app.command("compare")
def compare_command(
    reference: Annotated[
        Path | None,
        typer.Option("--reference", help="WT or reference analysis table."),
    ] = None,
    variant: Annotated[
        Path | None,
        typer.Option("--variant", help="PTM variant analysis table."),
    ] = None,
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Comparison output path."),
    ] = Path("comparison.csv"),
) -> None:
    """Placeholder for WT-vs-PTM comparisons."""

    typer.echo(
        f"Stage 1 placeholder: would compare reference={reference} and variant={variant} "
        f"into {output}."
    )


@app.command("report")
def report_command(
    config: ConfigOption = None,
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Report output path."),
    ] = Path("report.html"),
) -> None:
    """Placeholder for report and figure generation."""

    typer.echo(f"Stage 1 placeholder: would generate report {output} using config {config}.")


def main() -> None:
    """Console-script entry point."""

    app()


if __name__ == "__main__":
    main()
