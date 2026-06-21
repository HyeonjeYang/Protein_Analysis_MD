"""Report generation for analyzed WT-vs-PTM projects."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from idrptm.analysis.compare import ProjectComparison, compare_project, load_project_replicates
from idrptm.plotting.plots import (
    plot_delta_matrix,
    plot_distribution,
    plot_lines,
    plot_matrix,
    plot_ptm_sites,
    plot_summary_table,
    save_figure,
)


@dataclass(frozen=True)
class ReportPlan:
    """A planned report artifact."""

    output: Path
    title: str = "protein_analysis_md report"


@dataclass(frozen=True)
class ReportResult:
    """Generated report paths."""

    output_dir: Path
    report_path: Path
    figure_paths: tuple[Path, ...]
    comparison: ProjectComparison


def build_report_plan(output: str | Path, title: str = "protein_analysis_md report") -> ReportPlan:
    """Create a report plan."""

    return ReportPlan(output=Path(output), title=title)


def generate_report(project_dir: str | Path, output_dir: str | Path | None = None) -> ReportResult:
    """Generate figures and a Markdown report for a project directory."""

    project_path = Path(project_dir)
    root = Path(output_dir) if output_dir is not None else project_path / "report"
    figures = root / "figures"
    root.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    comparison = compare_project(project_path, output_dir=project_path / "comparison")
    runs = load_project_replicates(project_path)
    manifest = pd.DataFrame(_read_manifest(project_path / "manifest.csv"))
    smoothing = _report_smoothing(project_path)
    figure_paths: list[Path] = []

    rg_distribution = _distribution_table(runs, "rg", "rg")
    ree_distribution = _distribution_table(runs, "ree", "ree")
    figure_paths.extend(
        save_figure(
            plot_distribution(rg_distribution, "rg", "Rg (nm)", "Rg distribution"),
            figures / "rg_distribution",
        )
    )
    figure_paths.extend(
        save_figure(
            plot_distribution(ree_distribution, "ree", "Ree (nm)", "Ree distribution"),
            figures / "ree_distribution",
        )
    )

    figure_paths.extend(
        save_figure(
            plot_matrix(
                comparison.map_means[comparison.wt_condition],
                f"Contact map: {comparison.wt_condition}",
                "Contact probability (dimensionless)",
            ),
            figures / "contact_map_wt",
        )
    )
    for condition, delta_map in comparison.delta_maps.items():
        figure_paths.extend(
            save_figure(
                plot_matrix(
                    comparison.map_means[condition],
                    f"Contact map: {condition}",
                    "Contact probability (dimensionless)",
                ),
                figures / f"contact_map_{_slug(condition)}",
            )
        )
        figure_paths.extend(
            save_figure(
                plot_delta_matrix(delta_map, f"Delta contact map: {condition} - WT"),
                figures / f"delta_contact_map_{_slug(condition)}",
            )
        )

    ps_aggregate = pd.read_parquet(comparison.outputs["ps_aggregate_parquet"])
    scaling_aggregate = _scaling_aggregate(runs)
    show_raw_points = bool(smoothing.get("show_raw_points", True))
    show_smoothed_line = bool(smoothing.get("show_smoothed_line", True))
    ps_smooth = (
        "p_smooth_mean"
        if smoothing.get("use_smoothed_ps", True) and "p_smooth_mean" in ps_aggregate
        else None
    )
    rs_smooth = (
        "distance_smooth_mean"
        if smoothing.get("use_smoothed_rs", True)
        and "distance_smooth_mean" in scaling_aggregate
        else None
    )
    figure_paths.extend(
        save_figure(
            plot_lines(
                ps_aggregate,
                "s",
                "p_mean",
                "P(s) (dimensionless)",
                "Contact probability P(s)",
                smooth_y=ps_smooth,
                show_raw_points=show_raw_points,
                show_smoothed_line=show_smoothed_line,
            ),
            figures / "ps",
        )
    )
    figure_paths.extend(
        save_figure(
            plot_lines(
                scaling_aggregate,
                "s",
                "distance_mean",
                "R(s) (nm)",
                "Internal distance R(s)",
                smooth_y=rs_smooth,
                show_raw_points=show_raw_points,
                show_smoothed_line=show_smoothed_line,
            ),
            figures / "rs",
        )
    )
    figure_paths.extend(
        save_figure(plot_ptm_sites(manifest), figures / "ptm_sites")
    )
    figure_paths.extend(
        save_figure(plot_summary_table(comparison.summary_table), figures / "summary_table")
    )

    report_path = root / "report.md"
    report_path.write_text(
        _report_markdown(project_path, comparison, figure_paths, smoothing),
        encoding="utf-8",
    )
    return ReportResult(
        output_dir=root,
        report_path=report_path,
        figure_paths=tuple(figure_paths),
        comparison=comparison,
    )


def _distribution_table(
    runs: tuple[object, ...],
    dataframe_attr: str,
    value_column: str,
) -> pd.DataFrame:
    tables: list[pd.DataFrame] = []
    for run in runs:
        table = getattr(run, dataframe_attr)[[value_column]].copy()
        table["condition"] = run.condition
        table["replicate_id"] = run.replicate_id
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def _scaling_aggregate(runs: tuple[object, ...]) -> pd.DataFrame:
    tables: list[pd.DataFrame] = []
    for run in runs:
        raw_column = "mean_distance_nm" if "mean_distance_nm" in run.scaling else "distance"
        columns = ["s", raw_column]
        if "mean_distance_nm_smooth" in run.scaling:
            columns.append("mean_distance_nm_smooth")
        table = run.scaling[columns].copy()
        table = table.rename(columns={raw_column: "distance"})
        table["condition"] = run.condition
        table["replicate_id"] = run.replicate_id
        tables.append(table)
    combined = pd.concat(tables, ignore_index=True)
    rows = []
    for (condition, separation), group in combined.groupby(["condition", "s"], sort=True):
        row = {
            "condition": condition,
            "s": int(separation),
            "distance_mean": float(group["distance"].mean()),
        }
        if "mean_distance_nm_smooth" in group:
            smooth_values = group["mean_distance_nm_smooth"].dropna()
            if len(smooth_values):
                row["distance_smooth_mean"] = float(smooth_values.mean())
        rows.append(row)
    return pd.DataFrame(rows)


def _report_markdown(
    project_dir: Path,
    comparison: ProjectComparison,
    figure_paths: list[Path],
    smoothing: dict[str, object],
) -> str:
    lines = [
        "# protein_analysis_md report",
        "",
        f"Project directory: `{project_dir}`",
        f"WT condition: `{comparison.wt_condition}`",
        _smoothing_report_line(smoothing),
        "",
        "## Summary Table",
        "",
        _markdown_table(comparison.summary_table),
        "",
        "## Figures",
        "",
    ]
    for path in figure_paths:
        if path.suffix == ".png":
            relative = path.relative_to(path.parents[1])
            lines.extend(
                [
                    f"### {path.stem.replace('_', ' ').title()}",
                    "",
                    f"![{path.stem}]({relative})",
                    "",
                ]
            )
    lines.extend(
        [
            "## Comparison Outputs",
            "",
            f"- Summary CSV: `{comparison.outputs['summary_csv']}`",
            f"- Delta P(s): `{comparison.outputs['delta_ps_parquet']}`",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def _report_smoothing(project_dir: Path) -> dict[str, object]:
    resolved_path = project_dir / "config_resolved.json"
    if not resolved_path.is_file():
        return {}
    try:
        payload = json.loads(resolved_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    report = payload.get("report", {})
    if not isinstance(report, dict):
        return {}
    smoothing = report.get("smoothing", {})
    return dict(smoothing) if isinstance(smoothing, dict) else {}


def _smoothing_report_line(smoothing: dict[str, object]) -> str:
    if not smoothing or not smoothing.get("show_smoothing_metadata", False):
        return "Smoothing: raw curves are shown unless smoothed columns are explicitly available."
    ps = "smoothed" if smoothing.get("use_smoothed_ps", False) else "raw"
    rs = "smoothed" if smoothing.get("use_smoothed_rs", False) else "raw"
    return f"Smoothing: P(s) plot uses `{ps}` trend display; R(s) plot uses `{rs}` trend display."


def _read_manifest(manifest_path: Path) -> list[dict[str, str]]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _markdown_table(table: pd.DataFrame) -> str:
    columns = list(table.columns)
    rows = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in table.iterrows():
        rows.append("| " + " | ".join(_format_cell(row[column]) for column in columns) + " |")
    return "\n".join(rows)


def _format_cell(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.5g}"
    return str(value)


def _slug(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in value
    )
