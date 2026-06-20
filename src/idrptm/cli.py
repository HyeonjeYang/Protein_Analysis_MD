"""Command-line interface for idr-ptm-md."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, cast

import typer

from idrptm import __version__
from idrptm.analysis.compare import compare_project
from idrptm.analysis.pipeline import analyze_run_directory
from idrptm.calvados_adapter import prepare_from_config_file
from idrptm.design import design_from_config_file
from idrptm.plotting.report import generate_report

app = typer.Typer(
    name="idrptm",
    help="Prepare, run, and analyze CALVADOS-backed IDR/PTM MD workflows.",
    no_args_is_help=True,
)


ConfigOption = Annotated[
    Path | None,
    typer.Option("--config", "-c", help="Workflow YAML configuration file."),
]
TrajectoryReader = Literal["mdtraj", "mdanalysis"]


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
    run_dir: Annotated[
        Path,
        typer.Argument(
            help="Prepared CALVADOS run directory containing top.pdb and trajectory.dcd."
        ),
    ],
    config: ConfigOption = None,
    topology: Annotated[
        Path | None,
        typer.Option("--topology", help="Topology PDB path. Defaults to RUN_DIR/top.pdb."),
    ] = None,
    trajectory: Annotated[
        Path | None,
        typer.Option(
            "--trajectory",
            "-t",
            help="Trajectory DCD path. Defaults to RUN_DIR/trajectory.dcd.",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output-dir",
            "-o",
            help="Analysis output directory. Defaults to RUN_DIR/analysis.",
        ),
    ] = None,
    trajectory_reader: Annotated[
        str,
        typer.Option(
            "--trajectory-reader",
            help="Trajectory reader backend: mdtraj or mdanalysis.",
        ),
    ] = "mdtraj",
) -> None:
    """Analyze a prepared CALVADOS trajectory."""

    try:
        if trajectory_reader not in {"mdtraj", "mdanalysis"}:
            raise ValueError("Trajectory reader must be 'mdtraj' or 'mdanalysis'.")
        result = analyze_run_directory(
            run_dir,
            config_path=config,
            topology=topology,
            trajectory=trajectory,
            trajectory_reader=cast(TrajectoryReader, trajectory_reader),
            output_dir=output_dir,
        )
    except Exception as exc:
        typer.echo(f"Analyze failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote analysis outputs to: {result.output_dir}")
    typer.echo(f"Summary: {result.summary_json}")


@app.command("compare")
def compare_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Project directory containing manifest.csv and analyzed runs."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Comparison output directory."),
    ] = None,
) -> None:
    """Compare each PTM condition against the detected WT condition."""

    try:
        result = compare_project(project_dir, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Compare failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"WT condition: {result.wt_condition}")
    typer.echo(f"Wrote comparison outputs to: {result.output_dir}")
    typer.echo(f"Summary: {result.outputs['summary_csv']}")


@app.command("report")
def report_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Project directory containing manifest.csv and analyzed runs."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Report output directory."),
    ] = None,
) -> None:
    """Generate WT-vs-PTM comparison figures and a Markdown report."""

    try:
        result = generate_report(project_dir, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Report failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote report: {result.report_path}")
    typer.echo(f"Wrote {len(result.figure_paths)} figure file(s).")


def main() -> None:
    """Console-script entry point."""

    app()


if __name__ == "__main__":
    main()
