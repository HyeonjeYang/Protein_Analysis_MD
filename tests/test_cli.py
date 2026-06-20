from __future__ import annotations

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
