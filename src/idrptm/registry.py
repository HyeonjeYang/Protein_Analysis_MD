"""Optional SQLite project registry."""

from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path


def update_project_registry(project_dir: str | Path) -> Path:
    """Build or refresh ``project.db`` from manifest and run status files."""

    root = Path(project_dir)
    db_path = root / "project.db"
    manifest = root / "manifest.csv"
    with sqlite3.connect(db_path) as connection:
        connection.execute(
            "create table if not exists runs "
            "(run_id text primary key, ptm_state text, status text, run_dir text)"
        )
        connection.execute("delete from runs")
        if manifest.exists():
            with manifest.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    run_id = row.get("variant_id") or row.get("run_id") or ""
                    run_dir = (
                        (root / row["metadata_path"]).parent
                        if row.get("metadata_path")
                        else root
                    )
                    status = _run_status(run_dir)
                    connection.execute(
                        "insert or replace into runs(run_id, ptm_state, status, run_dir) "
                        "values (?, ?, ?, ?)",
                        (run_id, row.get("ptm_state", ""), status, str(run_dir)),
                    )
        connection.commit()
    return db_path


def list_runs(project_dir: str | Path) -> list[dict[str, object]]:
    """Return run rows from the project registry."""

    db_path = update_project_registry(project_dir)
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        return [dict(row) for row in connection.execute("select * from runs order by run_id")]


def summarize_registry(project_dir: str | Path) -> dict[str, int]:
    """Return run counts by status."""

    db_path = update_project_registry(project_dir)
    with sqlite3.connect(db_path) as connection:
        rows = connection.execute("select status, count(*) from runs group by status")
        return {str(status): int(count) for status, count in rows}


def query_runs(project_dir: str | Path, expression: str) -> list[dict[str, object]]:
    """Run a simple pandas-style query over registry rows."""

    import pandas as pd

    rows = list_runs(project_dir)
    table = pd.DataFrame(rows)
    if table.empty:
        return []
    return table.query(expression).to_dict(orient="records")


def _run_status(run_dir: Path) -> str:
    status_path = run_dir / "run_status.json"
    if not status_path.exists():
        return "prepared" if (run_dir / "run.py").exists() else "unknown"
    payload = json.loads(status_path.read_text(encoding="utf-8"))
    return str(payload.get("status", "unknown"))
