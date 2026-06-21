from __future__ import annotations

from pathlib import Path

from protein_analysis_md.environment_check import (
    LOCAL_INTERPRETATION,
    REMOTE_INTERPRETATION,
    collect_environment_info,
    interpret_environment,
)
from protein_analysis_md.repo_check import run_repo_check


def test_interpret_environment_detects_remote_ssh() -> None:
    env = {"SSH_CONNECTION": "1.1.1.1 123 2.2.2.2 22"}

    assert interpret_environment(env) == REMOTE_INTERPRETATION


def test_interpret_environment_detects_local_without_markers() -> None:
    assert interpret_environment({}) == LOCAL_INTERPRETATION


def test_collect_environment_info_is_optional_dependency_safe(tmp_path) -> None:
    info = collect_environment_info(cwd=tmp_path, env={})

    assert info["cwd"] == str(tmp_path)
    assert "calvados" in info
    assert "openmm" in info
    assert info["interpretation"] == LOCAL_INTERPRETATION


def test_repo_check_writes_markdown(tmp_path) -> None:
    (tmp_path / "LICENSE").write_text("GPL-3.0-only\n", encoding="utf-8")
    (tmp_path / "NOTICE.md").write_text("# Notice\n", encoding="utf-8")
    (tmp_path / "THIRD_PARTY_LICENSES.md").write_text("# Third Party\n", encoding="utf-8")
    (tmp_path / "CITATION.cff").write_text("cff-version: 1.2.0\n", encoding="utf-8")
    (tmp_path / "README.md").write_text("# Test\n", encoding="utf-8")
    (tmp_path / ".gitignore").write_text("\n".join(["__pycache__/", "*.dcd"]) + "\n")

    result = run_repo_check(tmp_path, output_path=Path("REPO_CHECK.md"))

    assert result.output_path.exists()
    assert "not legal advice" in result.markdown
    assert result.findings["readme"]["line_count"] == 1
