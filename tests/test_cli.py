from __future__ import annotations

from textwrap import dedent

import pytest


def test_cli_help_works() -> None:
    pytest.importorskip("typer")
    from typer.testing import CliRunner

    from idrptm.cli import app

    runner = CliRunner()

    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "prepare" in result.output
    assert "analyze" in result.output
    assert "report" in result.output
    assert "search-uniprot" in result.output
    assert "estimate-size" in result.output
    assert "compile" in result.output
    assert "wizard" in result.output
    assert "status" in result.output
    assert "env-check" in result.output
    assert "repo-check" in result.output
    assert "launch-local" in result.output


def test_pamd_wrapper_help_works() -> None:
    pytest.importorskip("typer")
    from typer.testing import CliRunner

    from protein_analysis_md.cli import app

    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "pamd" in result.output


def test_design_cli_generates_outputs(tmp_path) -> None:
    pytest.importorskip("typer")
    pytest.importorskip("pydantic")
    pytest.importorskip("yaml")
    from typer.testing import CliRunner

    from idrptm.cli import app

    config = tmp_path / "config.yaml"
    output_dir = tmp_path / "designs"
    config.write_text(
        dedent(
            """
            project: cli_design
            sequence:
              name: cli_seq
              sequence: "AST"
            ptm:
              mode: explicit
              include_wt: true
              sites:
                - position: 2
                  residue: S
                  ptm: pSer
            runner:
              work_dir: ignored_by_test
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(app, ["design", str(config), "--output-dir", str(output_dir)])

    assert result.exit_code == 0, result.output
    assert (output_dir / "manifest.csv").exists()
    assert (output_dir / "fasta" / "cli_seq__WT.fasta").exists()
    assert (output_dir / "fasta" / "cli_seq__pSer2.fasta").exists()
    assert (output_dir / "runs" / "cli_seq__pSer2" / "metadata.yaml").exists()


def test_prepare_cli_dry_run(tmp_path) -> None:
    pytest.importorskip("typer")
    pytest.importorskip("pydantic")
    pytest.importorskip("yaml")
    from typer.testing import CliRunner

    from idrptm.cli import app

    config = tmp_path / "config.yaml"
    output_dir = tmp_path / "prepare_dry_run"
    config.write_text(
        dedent(
            """
            project: cli_prepare
            sequence:
              name: prep_seq
              sequence: "AST"
            ptm:
              mode: wt
              include_wt: true
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        app,
        ["prepare", str(config), "--output-dir", str(output_dir), "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert "Would prepare" in result.output
    assert not output_dir.exists()


def test_analyze_cli_reports_missing_trajectory_files(tmp_path) -> None:
    pytest.importorskip("typer")
    from typer.testing import CliRunner

    from idrptm.cli import app

    result = CliRunner().invoke(app, ["analyze", str(tmp_path)])

    assert result.exit_code == 1
    assert "Missing topology PDB" in result.output
