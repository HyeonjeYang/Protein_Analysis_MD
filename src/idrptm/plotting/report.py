"""Report generation for analyzed WT-vs-PTM projects."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from idrptm.analysis.compare import ProjectComparison, compare_project, load_project_replicates
from idrptm.plotting.plots import (
    plot_contact_eigenvectors,
    plot_contact_loading_heatmap,
    plot_delta_ev,
    plot_delta_matrix,
    plot_distribution,
    plot_ev1_correlation,
    plot_explained_variance,
    plot_lines,
    plot_matrix,
    plot_nmf_residue_weights,
    plot_pca_centroid_shift,
    plot_pca_score_scatter,
    plot_pca_timeseries,
    plot_ptm_sites,
    plot_summary_table,
    save_figure,
)
from idrptm.visualization import VisualizationArtifact, save_visualization
from idrptm.visualization.cleavage import (
    event_schedule_figure,
    fragment_architecture_figure,
    plot_cut_number_trend,
    plot_fragment_length_distribution,
)
from idrptm.visualization.free_energy import plot_free_energy_grid
from idrptm.visualization.heatmaps import (
    plot_heatmap,
    ptm_site_contact_profile,
    residue_class_contact_matrix,
)
from idrptm.visualization.phase import (
    cluster_distribution_figure,
    phase_reliability_warning,
    plot_inter_chain_contact_heatmap,
)
from idrptm.visualization.ptm import (
    delta_ptm_site_profile,
    plot_ptm_site_profile,
    plot_residue_class_contact_changes,
    residue_class_contact_changes,
)
from idrptm.visualization.sequence_tracks import (
    plot_sequence_summary,
    plot_sequence_tracks,
    sequence_feature_table,
)
from idrptm.visualization.single_chain import (
    local_scaling_exponent,
    plot_contact_degree,
    plot_local_scaling_exponent,
    plot_observed_expected_contact_map,
    plot_rg_ree_hexbin,
    plot_rg_ree_timeseries,
    rg_ree_plotting_data,
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
    if "mean_square_distance_mean" in scaling_aggregate:
        figure_paths.extend(
            save_figure(
                plot_lines(
                    scaling_aggregate,
                    "s",
                    "mean_square_distance_mean",
                    "<R^2(s)> (nm^2)",
                    "Mean-square internal distance <R^2(s)>",
                    show_raw_points=False,
                    show_smoothed_line=False,
                ),
                figures / "r2s",
            )
        )
    figure_paths.extend(
        save_figure(plot_ptm_sites(manifest), figures / "ptm_sites")
    )
    figure_paths.extend(
        save_figure(plot_summary_table(comparison.summary_table), figures / "summary_table")
    )
    figure_paths.extend(_single_chain_figures(runs, manifest, figures))
    figure_paths.extend(_ptm_figures(comparison, manifest, figures))
    figure_paths.extend(_cleavage_figures(project_path, comparison.summary_table, figures))
    figure_paths.extend(_phase_figures(runs, figures))
    figure_paths.extend(_free_energy_figures(runs, figures))
    figure_paths.extend(_decomposition_figures(runs, figures, comparison))

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
        columns.extend(
            column
            for column in (
                "mean_distance_nm_smooth",
                "mean_square_distance_nm2",
                "rms_distance_nm",
            )
            if column in run.scaling
        )
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
        if "mean_square_distance_nm2" in group:
            square_values = group["mean_square_distance_nm2"].dropna()
            if len(square_values):
                row["mean_square_distance_mean"] = float(square_values.mean())
        if "rms_distance_nm" in group:
            rms_values = group["rms_distance_nm"].dropna()
            if len(rms_values):
                row["rms_distance_mean"] = float(rms_values.mean())
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
        (
            "Decomposition figures, when present, are exploratory contact-environment/PCA "
            "analyses and not chromosome compartments."
        ),
        (
            "Free-energy landscapes, when present, are sampling-dependent visual summaries; "
            "short trajectories should not be used to infer kinetic barriers."
        ),
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
            (
                "- Decomposition comparison: "
                f"`{comparison.outputs.get('decomposition_comparison_csv', 'not available')}`"
            ),
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


def _single_chain_figures(
    runs: tuple[object, ...],
    manifest: pd.DataFrame,
    figures: Path,
) -> list[Path]:
    first = runs[0]
    paths: list[Path] = []
    artifacts = [
        save_visualization(
            plot_rg_ree_timeseries(first.rg, first.ree),
            rg_ree_plotting_data(first.rg, first.ree),
            figures / "rg_ree_timeseries",
            metadata={"raw_data_preserved": True, "smoothing": "visual overlay only"},
        ),
        save_visualization(
            plot_rg_ree_hexbin(first.rg, first.ree),
            rg_ree_plotting_data(first.rg, first.ree),
            figures / "rg_ree_hexbin",
            metadata={"raw_data_preserved": True},
        ),
        save_visualization(
            plot_contact_degree(first.contact_map),
            pd.DataFrame(
                {
                    "residue_index": range(1, first.contact_map.shape[0] + 1),
                    "contact_degree": first.contact_map.sum(axis=1),
                }
            ),
            figures / "contact_degree",
            metadata={"raw_data_preserved": True},
        ),
        save_visualization(
            plot_observed_expected_contact_map(first.contact_map),
            first.contact_map,
            figures / "contact_observed_expected_wt",
            metadata={"raw_data_preserved": True, "expected_by_separation": True},
        ),
    ]
    for artifact in artifacts:
        paths.extend(_artifact_figures(artifact))
    local = local_scaling_exponent(first.scaling)
    if not local.empty:
        paths.extend(
            _artifact_figures(
                save_visualization(
                    plot_local_scaling_exponent(local),
                    local,
                    figures / "local_scaling_exponent",
                    metadata={"visualization_only": True},
                )
            )
        )
    sequence = str(manifest.iloc[0].get("original_sequence", ""))
    if sequence:
        paths.extend(
            _artifact_figures(
                save_visualization(
                    plot_sequence_tracks(sequence, ptm_sites=_parse_ptm_sites_manifest(manifest)),
                    sequence_feature_table(sequence),
                    figures / "sequence_tracks",
                    metadata={"raw_data_preserved": True},
                )
            )
        )
        paths.extend(
            _artifact_figures(
                save_visualization(
                    plot_sequence_summary(sequence),
                    sequence_feature_table(sequence),
                    figures / "sequence_summary",
                    metadata={"raw_data_preserved": True},
                )
            )
        )
    return paths


def _ptm_figures(
    comparison: ProjectComparison,
    manifest: pd.DataFrame,
    figures: Path,
) -> list[Path]:
    paths: list[Path] = []
    sequence = str(manifest.iloc[0].get("original_sequence", ""))
    for condition, contact_map in comparison.map_means.items():
        if condition == comparison.wt_condition:
            continue
        sites = _condition_ptm_sites(manifest, condition)
        if sites:
            profile = ptm_site_contact_profile(contact_map, sites)
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_ptm_site_profile(contact_map, sites, condition=condition),
                        profile,
                        figures / f"ptm_site_profile_{_slug(condition)}",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
            delta_profile = delta_ptm_site_profile(
                comparison.map_means[comparison.wt_condition],
                contact_map,
                sites,
            )
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_ptm_site_profile(
                            contact_map - comparison.map_means[comparison.wt_condition],
                            sites,
                            condition=f"{condition} - WT",
                        ),
                        delta_profile,
                        figures / f"delta_ptm_site_profile_{_slug(condition)}",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
        if sequence and len(sequence) == contact_map.shape[0]:
            class_delta = residue_class_contact_changes(
                comparison.map_means[comparison.wt_condition],
                contact_map,
                sequence,
            )
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_residue_class_contact_changes(class_delta),
                        class_delta,
                        figures / f"residue_class_contact_changes_{_slug(condition)}",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
            class_matrix = residue_class_contact_matrix(contact_map, sequence)
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_heatmap(
                            class_matrix,
                            title="Residue-class contact matrix",
                            colorbar_label="Contact probability (dimensionless)",
                            x_label="Residue class",
                            y_label="Residue class",
                        ),
                        class_matrix,
                        figures / f"residue_class_contact_matrix_{_slug(condition)}",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
    return paths


def _cleavage_figures(project_dir: Path, summary: pd.DataFrame, figures: Path) -> list[Path]:
    paths: list[Path] = []
    cleavage_manifest = project_dir / "cleavage_manifest.csv"
    if cleavage_manifest.is_file():
        fragments = pd.read_csv(cleavage_manifest)
        if {"original_start", "original_end", "fragment_id"}.issubset(fragments.columns):
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        fragment_architecture_figure(fragments),
                        fragments,
                        figures / "fragment_architecture",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_fragment_length_distribution(fragments),
                        fragments,
                        figures / "fragment_length_distribution",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
    cleavage_sites = project_dir / "cleavage_sites.csv"
    if cleavage_sites.is_file():
        sites = pd.read_csv(cleavage_sites)
        if {"event_time_ns", "cut_number"}.issubset(sites.columns):
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        event_schedule_figure(sites),
                        sites,
                        figures / "cleavage_event_schedule",
                        metadata={"raw_data_preserved": True, "smoothing": "disabled"},
                    )
                )
            )
    if "cut_number" in summary and "delta_mean_Rg" in summary:
        paths.extend(
            _artifact_figures(
                save_visualization(
                    plot_cut_number_trend(
                        summary.dropna(subset=["cut_number"]),
                        "delta_mean_Rg",
                        ylabel="Delta mean Rg (nm)",
                        title="Cut number vs Rg",
                    ),
                    summary,
                    figures / "cut_number_delta_rg",
                    metadata={"raw_points": True, "quasi_dynamic": True},
                )
            )
        )
    return paths


def _phase_figures(runs: tuple[object, ...], figures: Path) -> list[Path]:
    paths: list[Path] = []
    for run in runs:
        cluster_path = run.analysis_dir / "cluster_size.parquet"
        if cluster_path.is_file():
            cluster = pd.read_parquet(cluster_path)
            warning = phase_reliability_warning(
                len(cluster),
                int(cluster["largest_cluster_size"].max()),
            )
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        cluster_distribution_figure(cluster),
                        cluster,
                        figures / f"cluster_size_distribution_{_slug(run.condition)}",
                        metadata={"warning": warning} if warning else None,
                    )
                )
            )
        inter_chain_path = run.analysis_dir / "inter_chain_contact_map.npy"
        if inter_chain_path.is_file():
            import numpy as np

            matrix = np.load(inter_chain_path)
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_inter_chain_contact_heatmap(matrix),
                        matrix,
                        figures / f"inter_chain_contact_map_{_slug(run.condition)}",
                        metadata={"raw_data_preserved": True},
                    )
                )
            )
    return paths


def _free_energy_figures(runs: tuple[object, ...], figures: Path) -> list[Path]:
    paths: list[Path] = []
    for run in runs:
        for free_energy_path in sorted(run.analysis_dir.glob("free_energy_*_free_energy.npy")):
            import numpy as np

            prefix = free_energy_path.name.removesuffix("_free_energy.npy")
            x_edges = np.load(run.analysis_dir / f"{prefix}_x_edges.npy")
            y_edges = np.load(run.analysis_dir / f"{prefix}_y_edges.npy")
            grid = np.load(free_energy_path)
            paths.extend(
                _artifact_figures(
                    save_visualization(
                        plot_free_energy_grid(
                            grid,
                            x_edges,
                            y_edges,
                            x_label="Rg" if "rg" in prefix.lower() else "x",
                            y_label="Ree" if "ree" in prefix.lower() else "y",
                            title=prefix.replace("_", " "),
                        ),
                        grid,
                        figures / f"{prefix}_{_slug(run.condition)}",
                        metadata={"exploratory": True, "raw_counts_preserved": True},
                    )
                )
            )
    return paths


def _decomposition_figures(
    runs: tuple[object, ...],
    figures: Path,
    comparison: ProjectComparison,
) -> list[Path]:
    paths: list[Path] = []
    delta_ev_path = comparison.outputs.get("delta_ev_csv")
    if delta_ev_path and delta_ev_path.is_file():
        paths.extend(
            save_figure(
                plot_delta_ev(pd.read_csv(delta_ev_path)),
                figures / "delta_contact_environment_ev1",
            )
        )
    decomposition_path = comparison.outputs.get("decomposition_comparison_csv")
    if decomposition_path and decomposition_path.is_file():
        paths.extend(
            save_figure(
                plot_ev1_correlation(pd.read_csv(decomposition_path)),
                figures / "ev1_correlation_series",
            )
        )
    centroid_path = comparison.outputs.get("pca_centroid_shift_csv")
    if centroid_path and centroid_path.is_file():
        paths.extend(
            save_figure(
                plot_pca_centroid_shift(pd.read_csv(centroid_path)),
                figures / "pca_centroid_shift",
            )
        )
    first = next((run for run in runs if (run.analysis_dir / "contact_eigs.csv").is_file()), None)
    if first is not None:
        eigs = pd.read_csv(first.analysis_dir / "contact_eigs.csv")
        paths.extend(
            save_figure(
                plot_contact_eigenvectors(eigs),
                figures / "contact_environment_eigenvectors",
            )
        )
        oe_path = first.analysis_dir / "contact_oe.npy"
        corr_path = first.analysis_dir / "contact_correlation.npy"
        if oe_path.is_file():
            import numpy as np

            paths.extend(
                save_figure(
                    plot_matrix(
                        np.load(oe_path),
                        "Observed/expected-like contact map",
                        "log contact enrichment (a.u.)",
                    ),
                    figures / "contact_observed_expected",
                )
            )
        if corr_path.is_file():
            import numpy as np

            paths.extend(
                save_figure(
                    plot_matrix(
                        np.load(corr_path),
                        "Contact-environment correlation",
                        "Correlation (dimensionless)",
                    ),
                    figures / "contact_environment_correlation",
                )
            )
    for source in ("feature_pca", "contact_pca", "pca"):
        run = next(
            (item for item in runs if (item.analysis_dir / f"{source}_scores.parquet").is_file()),
            None,
        )
        if run is None:
            continue
        scores = pd.read_parquet(run.analysis_dir / f"{source}_scores.parquet")
        if {"PC1", "PC2"}.issubset(scores.columns):
            paths.extend(
                save_figure(
                    plot_pca_score_scatter(scores, title=f"{source} PC1/PC2 scores"),
                    figures / f"{source}_score_scatter",
                )
            )
            paths.extend(
                save_figure(
                    plot_pca_timeseries(scores, title=f"{source} score time series"),
                    figures / f"{source}_timeseries",
                )
            )
        variance = run.analysis_dir / f"{source}_explained_variance.csv"
        if variance.is_file():
            paths.extend(
                save_figure(
                    plot_explained_variance(pd.read_csv(variance), f"{source} explained variance"),
                    figures / f"{source}_explained_variance",
                )
            )
        if source == "contact_pca":
            loading_path = run.analysis_dir / "contact_pca_loadings.npy"
            top_contacts_path = run.analysis_dir / "top_loading_contacts.csv"
            if loading_path.is_file() and top_contacts_path.is_file():
                import numpy as np

                paths.extend(
                    save_figure(
                        plot_contact_loading_heatmap(
                            np.load(loading_path),
                            pd.read_csv(top_contacts_path),
                        ),
                        figures / "contact_pca_loading_heatmap",
                    )
                )
        nmf_weights = run.analysis_dir / "nmf_residue_weights.csv"
        if nmf_weights.is_file():
            paths.extend(
                save_figure(
                    plot_nmf_residue_weights(pd.read_csv(nmf_weights)),
                    figures / "nmf_residue_weights",
                )
            )
    return paths


def _artifact_figures(artifact: VisualizationArtifact) -> list[Path]:
    paths = [artifact.png]
    if artifact.pdf is not None:
        paths.append(artifact.pdf)
    return paths


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


def _parse_ptm_sites_manifest(manifest: pd.DataFrame) -> list[int]:
    sites: set[int] = set()
    for value in manifest.get("ptm_sites_1based", []):
        sites.update(_parse_ptm_sites(str(value)))
    return sorted(sites)


def _condition_ptm_sites(manifest: pd.DataFrame, condition: str) -> list[int]:
    sites: set[int] = set()
    condition_values = (
        manifest["condition"] if "condition" in manifest else pd.Series("", index=manifest.index)
    )
    ptm_values = (
        manifest["ptm_state"] if "ptm_state" in manifest else pd.Series("", index=manifest.index)
    )
    rows = manifest[(condition_values == condition) | (ptm_values == condition)]
    for value in rows.get("ptm_sites_1based", []):
        sites.update(_parse_ptm_sites(str(value)))
    return sorted(sites)


def _parse_ptm_sites(site_text: str) -> list[int]:
    positions: list[int] = []
    for match in site_text.split(";"):
        if not match:
            continue
        digits = "".join(character for character in match if character.isdigit())
        if digits:
            positions.append(int(digits))
    return positions


def _slug(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in value
    )
