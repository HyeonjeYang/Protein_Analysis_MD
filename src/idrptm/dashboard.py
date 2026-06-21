"""Generate a local static HTML dashboard for a project."""

from __future__ import annotations

import csv
import html
import json
import os
import time
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from idrptm.provenance import slugify
from idrptm.watch import ProjectWatchSnapshot, summarize_watch


@dataclass(frozen=True)
class DashboardResult:
    """Generated static dashboard paths."""

    output_dir: Path
    index_html: Path
    data_json: Path
    opened: bool = False


def generate_dashboard(
    project_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    title: str | None = None,
    refresh_seconds: int | None = None,
    open_browser: bool = False,
) -> DashboardResult:
    """Write a local browser-friendly dashboard for a project directory."""

    root = Path(project_dir).resolve()
    out = Path(output_dir).resolve() if output_dir is not None else root / "dashboard"
    out.mkdir(parents=True, exist_ok=True)
    snapshot = summarize_watch(root)
    payload = _dashboard_payload(root, snapshot)
    data_json = out / "dashboard_data.json"
    data_json.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    index = out / "index.html"
    index.write_text(
        _html_page(
            root=root,
            out=out,
            title=title or f"{root.name} dashboard",
            snapshot=snapshot,
            payload=payload,
            refresh_seconds=refresh_seconds,
        ),
        encoding="utf-8",
    )
    opened = False
    if open_browser:
        opened = webbrowser.open(index.resolve().as_uri())
    return DashboardResult(output_dir=out, index_html=index, data_json=data_json, opened=opened)


def _dashboard_payload(root: Path, snapshot: ProjectWatchSnapshot) -> dict[str, Any]:
    manifest_rows = _read_csv(root / "manifest.csv")
    run_rows = _run_rows(root, manifest_rows)
    figures = _files(root / "report" / "figures", suffixes=(".png", ".pdf"))
    comparison = _files(root / "comparison", suffixes=(".csv", ".json", ".parquet", ".npy"))
    pymol = _files(root / "pymol", suffixes=(".pml", ".csv", ".md"))
    bundles = sorted(root.glob("*_bundle.tar.gz"))
    return {
        "project_dir": str(root),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "watch": {
            "pipeline_status": snapshot.pipeline_status,
            "n_runs": snapshot.n_runs,
            "status_counts": snapshot.status_counts,
            "ready_trajectories": snapshot.ready_trajectories,
            "total_dcd_mb": snapshot.total_dcd_mb,
            "frame_progress": snapshot.frame_progress,
        },
        "manifest_rows": len(manifest_rows),
        "runs": run_rows,
        "figures": [str(path) for path in figures],
        "comparison": [str(path) for path in comparison],
        "pymol": [str(path) for path in pymol],
        "bundles": [str(path) for path in bundles],
    }


def _html_page(
    *,
    root: Path,
    out: Path,
    title: str,
    snapshot: ProjectWatchSnapshot,
    payload: dict[str, Any],
    refresh_seconds: int | None,
) -> str:
    done, total = snapshot.frame_progress or (0, 0)
    percent = (100.0 * done / total) if total else 0.0
    refresh = (
        f'<meta http-equiv="refresh" content="{int(refresh_seconds)}">\n'
        if refresh_seconds
        else ""
    )
    sections = [
        _summary_section(root, out, snapshot, percent),
        _action_section(root, out),
        _figures_section(root, out, [Path(path) for path in payload["figures"]]),
        _links_section(
            "Comparison Data",
            root,
            out,
            [Path(path) for path in payload["comparison"]],
        ),
        _links_section("PyMOL", root, out, [Path(path) for path in payload["pymol"]]),
        _links_section("Bundles", root, out, [Path(path) for path in payload["bundles"]]),
        _runs_section(root, out, payload["runs"]),
    ]
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
{refresh}<title>{_esc(title)}</title>
<style>
:root {{
  color-scheme: light;
  --ink: #17202a;
  --muted: #5d6975;
  --line: #d9e1e8;
  --surface: #ffffff;
  --band: #f6f8fa;
  --accent: #0b6bcb;
  --ok: #1b7f47;
  --warn: #9a6700;
  --bad: #b42318;
}}
body {{
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  color: var(--ink);
  background: var(--band);
}}
header {{
  padding: 24px 28px 18px;
  background: #12263a;
  color: white;
}}
h1 {{ margin: 0 0 6px; font-size: 28px; letter-spacing: 0; }}
h2 {{ margin: 0 0 12px; font-size: 18px; letter-spacing: 0; }}
a {{ color: var(--accent); text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
main {{ padding: 22px 28px 40px; display: grid; gap: 18px; }}
section {{
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 18px;
}}
.muted {{ color: var(--muted); }}
.metrics {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 10px;
}}
.metric {{
  border: 1px solid var(--line);
  border-radius: 6px;
  padding: 12px;
  background: #fbfcfd;
}}
.metric strong {{ display: block; font-size: 22px; margin-top: 4px; }}
.bar {{ height: 12px; background: #e7edf3; border-radius: 999px; overflow: hidden; }}
.bar span {{ display: block; height: 100%; background: var(--ok); }}
.commands code {{
  display: block;
  white-space: pre-wrap;
  background: #f0f3f6;
  border-radius: 6px;
  padding: 10px;
  margin: 8px 0;
}}
table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--line); text-align: left; }}
th {{ color: var(--muted); font-weight: 600; }}
.gallery {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 14px;
}}
.figure img {{ width: 100%; border: 1px solid var(--line); border-radius: 6px; background: white; }}
.pill {{
  display: inline-block;
  padding: 2px 7px;
  border-radius: 999px;
  background: #edf4ff;
  color: #174ea6;
}}
</style>
</head>
<body>
<header>
  <h1>{_esc(title)}</h1>
  <div class="muted">{_esc(str(root))} · generated {_esc(str(payload["generated_at"]))}</div>
</header>
<main>
{''.join(sections)}
</main>
</body>
</html>
"""


def _summary_section(root: Path, out: Path, snapshot: ProjectWatchSnapshot, percent: float) -> str:
    pipeline = snapshot.pipeline_status
    counts = ", ".join(
        f"{_esc(key)}={value}" for key, value in sorted(snapshot.status_counts.items())
    )
    done, total = snapshot.frame_progress or (0, 0)
    return f"""
<section>
  <h2>Project Status</h2>
  <div class="metrics">
    <div class="metric">Pipeline<strong>{_esc(str(pipeline.get("stage", "unknown")))}</strong></div>
    <div class="metric">State<strong>{_esc(str(pipeline.get("status", "unknown")))}</strong></div>
    <div class="metric">Runs<strong>{snapshot.n_runs}</strong></div>
    <div class="metric">
      DCD Ready<strong>{snapshot.ready_trajectories}/{snapshot.n_runs}</strong>
    </div>
    <div class="metric">DCD Size<strong>{snapshot.total_dcd_mb:.3f} MB</strong></div>
  </div>
  <p class="muted">Status counts: {counts or "none"}</p>
  <div class="bar"><span style="width:{percent:.1f}%"></span></div>
  <p class="muted">Frame estimate: {done}/{total} ({percent:.1f}%)</p>
  <p><a href="{_href(root / "pipeline_status.json", out)}">pipeline_status.json</a></p>
</section>
"""


def _action_section(root: Path, out: Path) -> str:
    return f"""
<section class="commands">
  <h2>Useful Commands</h2>
  <code>pamd watch {_esc(str(root))} --follow</code>
  <code>pamd finalize {_esc(str(root))}</code>
  <code>pamd dashboard {_esc(str(root))} --open</code>
  <code>pamd pack {_esc(str(root))}</code>
  <p><a href="{_href(out / "dashboard_data.json", out)}">dashboard_data.json</a></p>
</section>
"""


def _figures_section(root: Path, out: Path, figures: list[Path]) -> str:
    pngs = [path for path in figures if path.suffix.lower() == ".png"]
    if not pngs:
        return """
<section>
  <h2>Figures</h2>
  <p class="muted">
    No PNG figures found yet. Run <code>pamd finalize</code> after simulations finish.
  </p>
</section>
"""
    items = []
    for path in pngs[:24]:
        items.append(
            f"""
    <div class="figure">
      <a href="{_href(path, out)}"><img src="{_href(path, out)}" alt="{_esc(path.name)}"></a>
      <p>{_esc(path.name)}</p>
    </div>
"""
        )
    return f"""
<section>
  <h2>Figures</h2>
  <div class="gallery">{''.join(items)}</div>
  <p class="muted">Showing up to 24 PNG figures.
  Full figure folder: <a href="{_href(root / "report" / "figures", out)}">report/figures</a></p>
</section>
"""


def _links_section(title: str, root: Path, out: Path, paths: list[Path]) -> str:
    if not paths:
        return f"""
<section>
  <h2>{_esc(title)}</h2>
  <p class="muted">No files found yet.</p>
</section>
"""
    rows = "\n".join(
        (
            f'<tr><td><a href="{_href(path, out)}">'
            f"{_esc(path.relative_to(root).as_posix())}</a></td></tr>"
        )
        for path in paths[:80]
    )
    return f"""
<section>
  <h2>{_esc(title)}</h2>
  <table><tbody>{rows}</tbody></table>
</section>
"""


def _runs_section(root: Path, out: Path, runs: list[dict[str, str]]) -> str:
    if not runs:
        return """
<section>
  <h2>Runs</h2>
  <p class="muted">No runs found.</p>
</section>
"""
    rows = []
    for run in runs[:200]:
        rows.append(
            "<tr>"
            f"<td>{_esc(run['run_id'])}</td>"
            f"<td><span class=\"pill\">{_esc(run['status'])}</span></td>"
            f"<td>{_link(run.get('analysis'), out, 'analysis')}</td>"
            f"<td>{_link(run.get('trajectory'), out, 'DCD')}</td>"
            f"<td>{_link(run.get('pymol'), out, 'PyMOL')}</td>"
            f"<td>{_link(run.get('parameters'), out, 'parameters')}</td>"
            "</tr>"
        )
    return f"""
<section>
  <h2>Runs</h2>
  <table>
    <thead>
      <tr>
        <th>Run</th><th>Status</th><th>Analysis</th>
        <th>Trajectory</th><th>PyMOL</th><th>Parameters</th>
      </tr>
    </thead>
    <tbody>{''.join(rows)}</tbody>
  </table>
</section>
"""


def _run_rows(root: Path, manifest_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for row in manifest_rows:
        run_id = row.get("variant_id") or row.get("run_id") or "run"
        run_dir = _run_dir(root, row)
        status = _read_json(run_dir / "run_status.json").get("status", "prepared")
        pymol_dir = root / "pymol" / _slug_for_pymol(run_id)
        rows.append(
            {
                "run_id": run_id,
                "status": str(status),
                "analysis": str(run_dir / "analysis" / "summary.json")
                if (run_dir / "analysis" / "summary.json").exists()
                else "",
                "trajectory": str(_latest(run_dir, "*.dcd") or ""),
                "pymol": str(pymol_dir / "load.pml") if (pymol_dir / "load.pml").exists() else "",
                "parameters": str(run_dir / "parameters.txt")
                if (run_dir / "parameters.txt").exists()
                else "",
            }
        )
    return rows


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _run_dir(root: Path, row: dict[str, str]) -> Path:
    if row.get("metadata_path"):
        return (root / row["metadata_path"]).parent
    return root / "runs" / (row.get("variant_id") or row.get("run_id") or "")


def _files(root: Path, *, suffixes: tuple[str, ...]) -> list[Path]:
    if not root.exists():
        return []
    normalized = tuple(suffix.lower() for suffix in suffixes)
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in normalized
    )


def _latest(root: Path, pattern: str) -> Path | None:
    paths = sorted(root.glob(pattern), key=lambda path: path.stat().st_mtime)
    return paths[-1] if paths else None


def _link(path_text: str | None, out: Path, label: str) -> str:
    if not path_text:
        return ""
    path = Path(path_text)
    return f'<a href="{_href(path, out)}">{_esc(label)}</a>'


def _href(path: Path, out: Path) -> str:
    if not path.exists():
        return "#"
    try:
        rel = os.path.relpath(path.resolve(), out.resolve())
    except OSError:
        rel = str(path)
    return html.escape(Path(rel).as_posix(), quote=True)


def _slug_for_pymol(value: str) -> str:
    return slugify(value)


def _esc(value: object) -> str:
    return html.escape(str(value), quote=True)
