"""Report generation for analyzed WT-vs-PTM projects."""

from __future__ import annotations

import csv
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
    figure_paths.extend(
        save_figure(
            plot_lines(
                ps_aggregate,
                "s",
                "p_mean",
                "P(s) (dimensionless)",
                "Contact probability P(s)",
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
        _report_markdown(project_path, comparison, figure_paths),
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
        table = run.scaling[["s", "distance"]].copy()
        table["condition"] = run.condition
        table["replicate_id"] = run.replicate_id
        tables.append(table)
    combined = pd.concat(tables, ignore_index=True)
    rows = []
    for (condition, separation), group in combined.groupby(["condition", "s"], sort=True):
        rows.append(
            {
                "condition": condition,
                "s": int(separation),
                "distance_mean": float(group["distance"].mean()),
            }
        )
    return pd.DataFrame(rows)


def _report_markdown(
    project_dir: Path,
    comparison: ProjectComparison,
    figure_paths: list[Path],
) -> str:
    lines = [
        "# protein_analysis_md report",
        "",
        f"Project directory: `{project_dir}`",
        f"WT condition: `{comparison.wt_condition}`",
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
