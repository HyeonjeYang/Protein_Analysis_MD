"""Command-line interface for protein_analysis_md."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal, cast

import typer

from idrptm import __version__
from idrptm.runner import RunPhase

app = typer.Typer(
    name="pamd",
    help="Prepare, run, and analyze CALVADOS-backed protein/IDR MD workflows.",
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
        typer.Option("--version", help="Show the protein_analysis_md version and exit."),
    ] = False,
) -> None:
    """Top-level CLI callback."""

    if version:
        typer.echo(f"protein_analysis_md {__version__}")
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
        f"Would initialize a protein_analysis_md project at {output} "
        f"(force={force})."
    )


@app.command("search-uniprot")
def search_uniprot_command(
    query: Annotated[str, typer.Argument(help="UniProt search query.")],
    reviewed: Annotated[
        bool,
        typer.Option("--reviewed/--all", help="Restrict to Swiss-Prot reviewed entries."),
    ] = True,
    organism: Annotated[
        str | None,
        typer.Option("--organism", help="Organism scientific name."),
    ] = None,
    size: Annotated[int, typer.Option("--size", help="Number of candidates to show.")] = 10,
) -> None:
    """Search UniProt/Swiss-Prot and print ranked candidates."""

    from idrptm.uniprot import format_candidates, search_uniprot

    try:
        candidates = search_uniprot(query, reviewed=reviewed, organism=organism, size=size)
    except Exception as exc:
        typer.echo(f"UniProt search failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(format_candidates(candidates))


@app.command("fetch-sequence")
def fetch_sequence_command(
    query: Annotated[str, typer.Argument(help="UniProt query or accession.")],
    reviewed: Annotated[
        bool,
        typer.Option("--reviewed/--all", help="Restrict to Swiss-Prot reviewed entries."),
    ] = True,
    organism: Annotated[
        str | None,
        typer.Option("--organism", help="Organism scientific name."),
    ] = None,
    accession: Annotated[
        str | None,
        typer.Option("--accession", help="Explicit UniProt accession."),
    ] = None,
    interactive: Annotated[
        bool,
        typer.Option("--interactive", help="Prompt to select among ambiguous candidates."),
    ] = False,
    yes: Annotated[
        bool,
        typer.Option("--yes", "-y", help="Accept a strong unique candidate non-interactively."),
    ] = False,
    refresh: Annotated[
        bool,
        typer.Option("--refresh", help="Refresh cached sequence metadata."),
    ] = False,
    cache_dir: Annotated[
        Path,
        typer.Option("--cache-dir", help="Sequence cache directory."),
    ] = Path("data/sequences"),
) -> None:
    """Fetch one UniProt/Swiss-Prot sequence into the local cache."""

    from idrptm.uniprot import fetch_sequence

    try:
        result = fetch_sequence(
            query,
            accession=accession,
            reviewed=reviewed,
            organism=organism,
            interactive=interactive,
            yes=yes,
            refresh=refresh,
            cache_dir=cache_dir,
        )
    except Exception as exc:
        typer.echo(f"UniProt fetch failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Fetched {result.candidate.accession} {result.candidate.entry_name}")
    typer.echo(f"FASTA: {result.fasta_path}")
    typer.echo(f"Metadata: {result.metadata_path}")


@app.command("estimate-size")
def estimate_size_command(
    config: Annotated[Path, typer.Argument(help="Workflow YAML configuration file.")],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for storage_estimate.json."),
    ] = None,
) -> None:
    """Estimate trajectory storage for a workflow before running simulations."""

    from idrptm.storage import estimate_from_config_file, format_storage_table

    try:
        estimate = estimate_from_config_file(config, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Storage estimate failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(format_storage_table(estimate))


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

    from idrptm.design import design_from_config_file

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

    from idrptm.calvados_adapter import prepare_from_config_file

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
    target: Annotated[
        Path,
        typer.Argument(help="Prepared run directory, or project directory with --all-runs."),
    ] = Path("."),
    phase: Annotated[
        str,
        typer.Option("--phase", help="Run phase: equilibration, production, or all."),
    ] = "all",
    all_runs: Annotated[
        bool,
        typer.Option("--all-runs", help="Run every prepared run directory in a project."),
    ] = False,
    python_executable: Annotated[
        str,
        typer.Option("--python", help="Python executable used to invoke generated run scripts."),
    ] = "python",
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run/--execute",
            help="Print planned commands instead of executing them.",
        ),
    ] = True,
) -> None:
    """Run prepared CALVADOS directories locally through generated run scripts."""

    from idrptm.runner import execute_local_runs, plan_local_runs, write_planned_status

    try:
        if phase not in {"equilibration", "production", "all"}:
            raise ValueError("Phase must be 'equilibration', 'production', or 'all'.")
        plans = plan_local_runs(
            target,
            phase=cast(RunPhase, phase),
            all_runs=all_runs,
            python_executable=python_executable,
        )
        if dry_run:
            typer.echo(f"Dry run: planned {len(plans)} local run(s).")
            for plan in plans:
                write_planned_status(plan)
                typer.echo(f"{plan.run_dir}: {' '.join(plan.command)}")
            return
        results = execute_local_runs(plans)
    except Exception as exc:
        typer.echo(f"Run failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    failed = [result for result in results if result.status == "failed"]
    for result in results:
        typer.echo(f"{result.run_dir}: {result.status} ({result.status_json})")
    if failed:
        raise typer.Exit(1)


@app.command("hpc-script")
def hpc_script_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Prepared project directory containing manifest.csv."),
    ],
    scheduler: Annotated[
        str,
        typer.Option("--scheduler", help="Scheduler script type. Currently only slurm."),
    ] = "slurm",
    phase: Annotated[
        str,
        typer.Option("--phase", help="Run phase: equilibration, production, or all."),
    ] = "all",
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for generated scheduler files."),
    ] = None,
    python_executable: Annotated[
        str,
        typer.Option("--python", help="Python executable used by array tasks."),
    ] = "python",
    submit: Annotated[
        bool,
        typer.Option("--submit", help="Submit the generated script. Not implemented yet."),
    ] = False,
) -> None:
    """Generate an HPC scheduler script without submitting it."""

    from idrptm.hpc import write_slurm_array_script

    try:
        if scheduler != "slurm":
            raise ValueError("Only --scheduler slurm is supported.")
        if phase not in {"equilibration", "production", "all"}:
            raise ValueError("Phase must be 'equilibration', 'production', or 'all'.")
        if submit:
            raise ValueError("Submission is intentionally not automatic yet; run sbatch manually.")
        script = write_slurm_array_script(
            project_dir,
            output_dir=output_dir,
            phase=cast(RunPhase, phase),
            python_executable=python_executable,
        )
    except Exception as exc:
        typer.echo(f"HPC script generation failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Wrote SLURM array script: {script}")


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

    from idrptm.analysis.pipeline import analyze_run_directory

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

    from idrptm.analysis.compare import compare_project

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

    from idrptm.plotting.report import generate_report

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
