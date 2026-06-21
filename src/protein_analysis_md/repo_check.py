"""Public repository readiness checks."""

from __future__ import annotations

import importlib
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REQUIRED_GITIGNORE_PATTERNS = (
    "__pycache__/",
    "*.py[cod]",
    ".pytest_cache/",
    ".ruff_cache/",
    ".mypy_cache/",
    ".ipynb_checkpoints/",
    ".venv/",
    "*.env",
    "build/",
    "dist/",
    "*.egg-info/",
    "runs/",
    "outputs/",
    "trajectories/",
    "analysis_outputs/",
    "reports/",
    "checkpoints/",
    "logs/",
    "*.dcd",
    "*.xtc",
    "*.trr",
    "*.nc",
    "*.h5",
    "*.chk",
    "*.cpt",
    "*.log",
    "*.out",
    "*.err",
    "*.npy",
    "*.npz",
    "*.parquet",
    "*.pkl",
    "*.sqlite",
    "*.db",
)

UPSTREAM_PATTERNS = (
    "CALVADOS",
    "OpenMM",
    "KULL",
    "Langevin",
    "Config(",
    "Components",
    "Simulation",
    "GPL",
    "MIT License",
    "Apache",
    "Copyright",
    "Adapted from",
)

SECRET_PATTERNS = (
    re.compile(r"OPENAI_API_KEY\s*[:=]", re.IGNORECASE),
    re.compile(r"HF_TOKEN\s*[:=]", re.IGNORECASE),
    re.compile(r"AWS_SECRET\w*\s*[:=]", re.IGNORECASE),
    re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----", re.IGNORECASE),
    re.compile(r"password\s*=", re.IGNORECASE),
)

GENERATED_SUFFIXES = (
    ".dcd",
    ".xtc",
    ".trr",
    ".nc",
    ".h5",
    ".chk",
    ".cpt",
    ".log",
    ".out",
    ".err",
    ".npy",
    ".npz",
    ".parquet",
    ".pkl",
    ".sqlite",
    ".db",
)

SKIP_DIRS = {
    ".git",
    ".venv",
    "venv",
    "env",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    "htmlcov",
    "work/calvados_venv",
}


@dataclass(frozen=True)
class RepoCheckResult:
    """Repository check result and Markdown report."""

    root: Path
    output_path: Path
    markdown: str
    findings: dict[str, Any]


def run_repo_check(
    root: str | Path = ".",
    *,
    output_path: str | Path = "REPO_CHECK.md",
) -> RepoCheckResult:
    """Run public-release checks and write a Markdown report."""

    root_path = Path(root).resolve()
    findings = {
        "license": _file_status(root_path, "LICENSE"),
        "notice": _file_status(root_path, "NOTICE.md"),
        "third_party": _file_status(root_path, "THIRD_PARTY_LICENSES.md"),
        "citation": _file_status(root_path, "CITATION.cff"),
        "readme": _readme_status(root_path),
        "gitignore": _gitignore_status(root_path),
        "large_files": _large_files(root_path),
        "generated_files": _generated_files(root_path),
        "upstream_references": _text_hits(root_path, UPSTREAM_PATTERNS),
        "secrets": _secret_hits(root_path),
        "imports": _import_status(),
        "tracked_generated": _tracked_generated(root_path),
    }
    markdown = _markdown(findings)
    output = root_path / output_path
    output.write_text(markdown, encoding="utf-8")
    return RepoCheckResult(root=root_path, output_path=output, markdown=markdown, findings=findings)


def _file_status(root: Path, name: str) -> dict[str, Any]:
    path = root / name
    return {
        "path": name,
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }


def _readme_status(root: Path) -> dict[str, Any]:
    path = root / "README.md"
    if not path.is_file():
        return {"exists": False, "line_count": 0, "within_120_lines": False}
    line_count = len(path.read_text(encoding="utf-8").splitlines())
    return {"exists": True, "line_count": line_count, "within_120_lines": line_count <= 120}


def _gitignore_status(root: Path) -> dict[str, Any]:
    path = root / ".gitignore"
    if not path.is_file():
        return {"exists": False, "missing_patterns": list(REQUIRED_GITIGNORE_PATTERNS)}
    text = path.read_text(encoding="utf-8")
    missing = [pattern for pattern in REQUIRED_GITIGNORE_PATTERNS if pattern not in text]
    return {"exists": True, "missing_patterns": missing}


def _large_files(root: Path) -> dict[str, list[str]]:
    files_5mb: list[str] = []
    files_50mb: list[str] = []
    for path in _iter_files(root):
        try:
            size = path.stat().st_size
        except OSError:
            continue
        rel = str(path.relative_to(root))
        if size > 5 * 1024 * 1024:
            files_5mb.append(f"{rel} ({size} bytes)")
        if size > 50 * 1024 * 1024:
            files_50mb.append(f"{rel} ({size} bytes)")
    return {"over_5mb": files_5mb[:100], "over_50mb": files_50mb[:100]}


def _generated_files(root: Path) -> list[str]:
    hits: list[str] = []
    for path in _iter_files(root):
        rel = str(path.relative_to(root))
        if path.suffix in GENERATED_SUFFIXES or rel.startswith(("runs/", "outputs/", "reports/")):
            hits.append(rel)
    return hits[:200]


def _tracked_generated(root: Path) -> list[str]:
    try:
        completed = subprocess.run(
            ["git", "ls-files"],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return []
    if completed.returncode != 0:
        return []
    tracked = []
    for name in completed.stdout.splitlines():
        path = Path(name)
        if path.suffix in GENERATED_SUFFIXES or name.startswith(("runs/", "outputs/", "reports/")):
            tracked.append(name)
    return tracked


def _text_hits(root: Path, patterns: tuple[str, ...]) -> list[str]:
    hits: list[str] = []
    lowered = [(pattern, pattern.lower()) for pattern in patterns]
    for path in _iter_text_files(root):
        text = _safe_read(path)
        if text is None:
            continue
        lower = text.lower()
        for pattern, lowered_pattern in lowered:
            if lowered_pattern in lower:
                hits.append(f"{path.relative_to(root)}: {pattern}")
                break
    return hits[:200]


def _secret_hits(root: Path) -> list[str]:
    hits: list[str] = []
    for path in _iter_text_files(root):
        text = _safe_read(path)
        if text is None:
            continue
        for pattern in SECRET_PATTERNS:
            if pattern.search(text):
                hits.append(f"{path.relative_to(root)}: {pattern.pattern}")
                break
    return hits[:100]


def _import_status() -> dict[str, Any]:
    modules = ("protein_analysis_md", "idrptm")
    status: dict[str, Any] = {}
    for module in modules:
        try:
            imported = importlib.import_module(module)
        except Exception as exc:
            status[module] = {"ok": False, "error": str(exc)}
        else:
            status[module] = {"ok": True, "version": getattr(imported, "__version__", None)}
    return status


def _markdown(findings: dict[str, Any]) -> str:
    lines = [
        "# Repository Check",
        "",
        "This report is a practical readiness scan, not legal advice. It warns about",
        "possible public-release issues but does not claim legal certainty.",
        "",
        "## Required Files",
        "",
    ]
    for key in ("license", "notice", "third_party", "citation"):
        value = findings[key]
        lines.append(f"- {value['path']}: {'present' if value['exists'] else 'missing'}")
    readme = findings["readme"]
    lines.extend(
        [
            f"- README.md lines: {readme['line_count']}",
            f"- README <= 120 lines: {readme['within_120_lines']}",
            "",
            "## Gitignore",
            "",
        ]
    )
    missing = findings["gitignore"]["missing_patterns"]
    lines.append("- Missing required patterns: " + (", ".join(missing) if missing else "none"))
    lines.extend(
        [
            "",
            "## Large Files",
            "",
            _bullet_list("Files > 5 MB", findings["large_files"]["over_5mb"]),
            _bullet_list("Files > 50 MB", findings["large_files"]["over_50mb"]),
            "",
            "## Generated Or MD Output Files",
            "",
            _bullet_list("Generated files found in working tree", findings["generated_files"]),
            _bullet_list("Tracked generated files", findings["tracked_generated"]),
            "",
            "## Possible Upstream/License References",
            "",
            _bullet_list("Reference hits", findings["upstream_references"]),
            "",
            "## Secret Scan",
            "",
            _bullet_list("Potential secret-pattern hits", findings["secrets"]),
            "",
            "## Import Check",
            "",
        ]
    )
    for module, status in findings["imports"].items():
        lines.append(f"- {module}: {'ok' if status['ok'] else 'failed'}")
    lines.extend(
        [
            "",
            "## Suggested Next Command",
            "",
            "```bash",
            "pytest",
            "ruff check .",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def _bullet_list(title: str, values: list[str]) -> str:
    if not values:
        return f"### {title}\n\n- none\n"
    rows = [f"### {title}", ""]
    rows.extend(f"- `{value}`" for value in values)
    rows.append("")
    return "\n".join(rows)


def _iter_files(root: Path):
    for path in root.rglob("*"):
        if path.is_file() and not _is_skipped(path, root):
            yield path


def _iter_text_files(root: Path):
    for path in _iter_files(root):
        if path.stat().st_size <= 2_000_000:
            yield path


def _safe_read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None


def _is_skipped(path: Path, root: Path) -> bool:
    rel = path.relative_to(root)
    parts = rel.parts
    for index in range(len(parts)):
        partial = "/".join(parts[: index + 1])
        if partial in SKIP_DIRS or parts[index] in SKIP_DIRS:
            return True
    return False
