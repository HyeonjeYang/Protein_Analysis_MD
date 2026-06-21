"""WT-vs-PTM comparison for analyzed project directories."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from idrptm.analysis.decomposition import ev1_correlation
from idrptm.analysis.scaling import fit_flory_exponent


@dataclass(frozen=True)
class RunAnalysis:
    """Loaded analysis outputs for one manifest row or replicate."""

    variant_id: str
    condition: str
    replicate_id: str
    ptm_sites_1based: str
    original_sequence: str
    cleavage_state: str
    cut_number: int | None
    event_time_ns: float | None
    analysis_dir: Path
    rg: pd.DataFrame
    ree: pd.DataFrame
    contact_map: np.ndarray
    ps: pd.DataFrame
    scaling: pd.DataFrame
    mean_rg: float
    mean_ree: float
    flory_exponent: float


@dataclass(frozen=True)
class ProjectComparison:
    """Paths and tables produced by a project comparison."""

    project_dir: Path
    output_dir: Path
    wt_condition: str
    summary_table: pd.DataFrame
    ps_delta: pd.DataFrame
    map_means: dict[str, np.ndarray]
    delta_maps: dict[str, np.ndarray]
    outputs: dict[str, Path]


@dataclass(frozen=True)
class ComparisonResult:
    """Small scalar comparison result retained for simple caller use."""

    reference: str
    variant: str
    metric: str
    delta: float | None = None


def compare_project(
    project_dir: str | Path,
    output_dir: str | Path | None = None,
) -> ProjectComparison:
    """Compare every PTM condition in ``project_dir`` against the detected WT."""

    project_path = Path(project_dir)
    runs = load_project_replicates(project_path)
    wt_condition = detect_wt_condition([run.condition for run in runs])
    root = Path(output_dir) if output_dir is not None else project_path / "comparison"
    root.mkdir(parents=True, exist_ok=True)

    summary = _summary_table(runs, wt_condition)
    ps_aggregate = _aggregate_ps(runs)
    ps_delta = _delta_ps(ps_aggregate, wt_condition)
    map_means = _aggregate_contact_maps(runs)
    delta_maps = {
        condition: matrix - map_means[wt_condition]
        for condition, matrix in map_means.items()
        if condition != wt_condition
    }

    outputs: dict[str, Path] = {}
    outputs["summary_csv"] = root / "comparison_summary.csv"
    outputs["summary_parquet"] = root / "comparison_summary.parquet"
    outputs["ps_aggregate_parquet"] = root / "ps_aggregate.parquet"
    outputs["delta_ps_parquet"] = root / "delta_ps.parquet"
    summary.to_csv(outputs["summary_csv"], index=False)
    summary.to_parquet(outputs["summary_parquet"], index=False)
    ps_aggregate.to_parquet(outputs["ps_aggregate_parquet"], index=False)
    ps_delta.to_parquet(outputs["delta_ps_parquet"], index=False)

    for condition, matrix in map_means.items():
        path = root / f"contact_map_mean_{_slug(condition)}.npy"
        np.save(path, matrix)
        outputs[f"contact_map_mean_{condition}"] = path
    for condition, matrix in delta_maps.items():
        path = root / f"delta_contact_map_{_slug(condition)}_minus_{_slug(wt_condition)}.npy"
        np.save(path, matrix)
        outputs[f"delta_contact_map_{condition}"] = path

    outputs.update(compare_decomposition_outputs(runs, wt_condition, root))

    metadata_path = root / "comparison_metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "project_dir": str(project_path),
                "wt_condition": wt_condition,
                "conditions": sorted({run.condition for run in runs}),
                "n_runs": len(runs),
                "outputs": {name: str(path) for name, path in outputs.items()},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["metadata_json"] = metadata_path

    return ProjectComparison(
        project_dir=project_path,
        output_dir=root,
        wt_condition=wt_condition,
        summary_table=summary,
        ps_delta=ps_delta,
        map_means=map_means,
        delta_maps=delta_maps,
        outputs=outputs,
    )


def load_project_replicates(project_dir: str | Path) -> tuple[RunAnalysis, ...]:
    """Load all analyzed manifest rows from a project directory."""

    project_path = Path(project_dir)
    manifest_path = project_path / "manifest.csv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing project manifest: {manifest_path}")

    rows = _read_manifest(manifest_path)
    analyses = tuple(_load_run_analysis(project_path, row) for row in rows)
    if not analyses:
        raise ValueError(f"No runs found in manifest: {manifest_path}")
    return analyses


def detect_wt_condition(conditions: list[str] | tuple[str, ...]) -> str:
    """Detect the WT condition from condition names."""

    unique_conditions = sorted(set(conditions))
    matches = [condition for condition in unique_conditions if is_wt_name(condition)]
    if not matches:
        raise ValueError("Could not detect a WT condition by name.")
    if len(matches) > 1:
        raise ValueError(f"Multiple WT-like conditions detected: {matches}")
    return matches[0]


def is_wt_name(name: str) -> bool:
    """Return true when a condition or variant name denotes wild type."""

    upper = name.upper()
    return upper == "WT" or re.search(r"(^|[_\-\s])WT($|[_\-\s])", upper) is not None


def compare_observable(reference: str, variant: str, metric: str) -> ComparisonResult:
    """Create a scalar comparison shell for callers that compare externally."""

    return ComparisonResult(reference=reference, variant=variant, metric=metric)


def compare_decomposition_outputs(
    runs: tuple[RunAnalysis, ...],
    wt_condition: str,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Compare optional decomposition outputs across conditions."""

    root = Path(output_dir)
    outputs: dict[str, Path] = {}
    ev_outputs = _compare_contact_eigs(runs, wt_condition)
    if ev_outputs:
        summary, delta = ev_outputs
        summary_path = root / "decomposition_comparison.csv"
        delta_path = root / "delta_ev.csv"
        summary.to_csv(summary_path, index=False)
        delta.to_csv(delta_path, index=False)
        outputs["decomposition_comparison_csv"] = summary_path
        outputs["delta_ev_csv"] = delta_path
    centroid = _compare_pca_centroids(runs, wt_condition)
    if not centroid.empty:
        path = root / "pca_centroid_shift.csv"
        centroid.to_csv(path, index=False)
        outputs["pca_centroid_shift_csv"] = path
    return outputs


def _read_manifest(manifest_path: Path) -> list[dict[str, str]]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_run_analysis(project_dir: Path, row: dict[str, str]) -> RunAnalysis:
    variant_id = row.get("variant_id") or row.get("name")
    if not variant_id:
        raise ValueError("Manifest row is missing variant_id.")
    condition = row.get("condition") or row.get("ptm_state") or variant_id
    replicate_id = row.get("replicate_id") or row.get("replicate") or variant_id
    cleavage_metadata = _cleavage_metadata(row)
    analysis_dir = _analysis_dir(project_dir, row, variant_id)
    _require_analysis_files(analysis_dir)

    rg = pd.read_parquet(analysis_dir / "timeseries_rg.parquet")
    ree = pd.read_parquet(analysis_dir / "timeseries_ree.parquet")
    contact_map = np.load(analysis_dir / "contact_map.npy")
    ps = pd.read_parquet(analysis_dir / "ps.parquet")
    scaling = pd.read_parquet(analysis_dir / "scaling.parquet")
    summary = _read_summary(analysis_dir / "summary.json")
    exponent = _flory_exponent(summary, scaling)

    return RunAnalysis(
        variant_id=variant_id,
        condition=condition,
        replicate_id=replicate_id,
        ptm_sites_1based=row.get("ptm_sites_1based", ""),
        original_sequence=row.get("original_sequence", ""),
        cleavage_state=str(cleavage_metadata["cleavage_state"]),
        cut_number=cleavage_metadata["cut_number"],
        event_time_ns=cleavage_metadata["event_time_ns"],
        analysis_dir=analysis_dir,
        rg=rg,
        ree=ree,
        contact_map=contact_map,
        ps=ps,
        scaling=scaling,
        mean_rg=float(rg["rg"].mean()),
        mean_ree=float(ree["ree"].mean()),
        flory_exponent=exponent,
    )


def _cleavage_metadata(row: dict[str, str]) -> dict[str, str | int | float | None]:
    components = _components_from_row(row)
    states = sorted(
        {
            str(component.get("cleavage_state"))
            for component in components
            if component.get("cleavage_state")
        }
    )
    cut_numbers = [
        int(component["cut_number"])
        for component in components
        if component.get("cut_number") not in {None, ""}
    ]
    event_times = [
        float(component["event_time_ns"])
        for component in components
        if component.get("event_time_ns") not in {None, ""}
    ]
    return {
        "cleavage_state": ";".join(states),
        "cut_number": max(cut_numbers) if cut_numbers else None,
        "event_time_ns": max(event_times) if event_times else None,
    }


def _components_from_row(row: dict[str, str]) -> list[dict[str, object]]:
    payload = row.get("components_json")
    if not payload:
        return []
    try:
        parsed = json.loads(payload)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    return [item for item in parsed if isinstance(item, dict)]


def _analysis_dir(project_dir: Path, row: dict[str, str], variant_id: str) -> Path:
    if row.get("analysis_path"):
        return project_dir / row["analysis_path"]
    if row.get("metadata_path"):
        return (project_dir / row["metadata_path"]).parent / "analysis"
    return project_dir / "runs" / variant_id / "analysis"


def _require_analysis_files(analysis_dir: Path) -> None:
    required = [
        "timeseries_rg.parquet",
        "timeseries_ree.parquet",
        "contact_map.npy",
        "ps.parquet",
        "scaling.parquet",
        "summary.json",
    ]
    missing = [name for name in required if not (analysis_dir / name).is_file()]
    if missing:
        missing_text = ", ".join(missing)
        raise FileNotFoundError(f"Missing analysis output(s) in {analysis_dir}: {missing_text}")


def _read_summary(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _flory_exponent(summary: dict[str, object], scaling: pd.DataFrame) -> float:
    fit = summary.get("flory_fit")
    if isinstance(fit, dict) and fit.get("nu") is not None:
        return float(fit["nu"])
    return fit_flory_exponent(scaling).nu


def _summary_table(runs: tuple[RunAnalysis, ...], wt_condition: str) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    by_condition = _group_runs(runs)
    wt_stats = _condition_scalar_stats(by_condition[wt_condition])
    for condition in sorted(by_condition):
        stats = _condition_scalar_stats(by_condition[condition])
        rows.append(
            {
                "condition": condition,
                "wt_condition": wt_condition,
                "n_replicates": len(by_condition[condition]),
                "mean_Rg": stats["mean_rg_mean"],
                "std_Rg": stats["mean_rg_std"],
                "sem_Rg": stats["mean_rg_sem"],
                "delta_mean_Rg": stats["mean_rg_mean"] - wt_stats["mean_rg_mean"],
                "mean_Ree": stats["mean_ree_mean"],
                "std_Ree": stats["mean_ree_std"],
                "sem_Ree": stats["mean_ree_sem"],
                "delta_mean_Ree": stats["mean_ree_mean"] - wt_stats["mean_ree_mean"],
                "flory_exponent": stats["flory_mean"],
                "std_flory_exponent": stats["flory_std"],
                "sem_flory_exponent": stats["flory_sem"],
                "delta_flory_exponent": stats["flory_mean"] - wt_stats["flory_mean"],
            }
        )
    return pd.DataFrame(rows)


def _condition_scalar_stats(runs: tuple[RunAnalysis, ...]) -> dict[str, float]:
    values = {
        "mean_rg": np.array([run.mean_rg for run in runs], dtype=float),
        "mean_ree": np.array([run.mean_ree for run in runs], dtype=float),
        "flory": np.array([run.flory_exponent for run in runs], dtype=float),
    }
    stats: dict[str, float] = {}
    for name, array in values.items():
        stats[f"{name}_mean"] = float(array.mean())
        stats[f"{name}_std"] = float(array.std(ddof=1)) if array.size > 1 else 0.0
        stats[f"{name}_sem"] = stats[f"{name}_std"] / float(np.sqrt(array.size))
    return stats


def _aggregate_ps(runs: tuple[RunAnalysis, ...]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for run in runs:
        table = run.ps.copy()
        if "p_contact" not in table and "p" in table:
            table["p_contact"] = table["p"]
        table["condition"] = run.condition
        table["replicate_id"] = run.replicate_id
        rows.append(table)
    combined = pd.concat(rows, ignore_index=True)

    aggregate_rows: list[dict[str, float | int | str]] = []
    for (condition, separation), group in combined.groupby(["condition", "s"], sort=True):
        values = group["p_contact"].to_numpy(dtype=float)
        std = float(values.std(ddof=1)) if values.size > 1 else 0.0
        row: dict[str, float | int | str] = {
            "condition": condition,
            "s": int(separation),
            "p_mean": float(values.mean()),
            "p_std": std,
            "p_sem": std / float(np.sqrt(values.size)),
            "n_replicates": int(values.size),
        }
        if "p_contact_smooth" in group:
            smooth_values = group["p_contact_smooth"].to_numpy(dtype=float)
            smooth_values = smooth_values[np.isfinite(smooth_values)]
            if smooth_values.size:
                smooth_std = (
                    float(smooth_values.std(ddof=1)) if smooth_values.size > 1 else 0.0
                )
                row["p_smooth_mean"] = float(smooth_values.mean())
                row["p_smooth_std"] = smooth_std
                row["p_smooth_sem"] = smooth_std / float(np.sqrt(smooth_values.size))
        aggregate_rows.append(row)
    return pd.DataFrame(aggregate_rows)


def _delta_ps(ps_aggregate: pd.DataFrame, wt_condition: str) -> pd.DataFrame:
    wt = ps_aggregate[ps_aggregate["condition"] == wt_condition]
    wt = wt.rename(columns={"p_mean": "p_wt_mean", "p_sem": "p_wt_sem"})
    rows = ps_aggregate[ps_aggregate["condition"] != wt_condition]
    merged = rows.merge(wt[["s", "p_wt_mean", "p_wt_sem"]], on="s", how="left")
    merged["delta_p"] = merged["p_mean"] - merged["p_wt_mean"]
    merged["delta_p_sem"] = np.sqrt(merged["p_sem"] ** 2 + merged["p_wt_sem"] ** 2)
    if {"p_smooth_mean", "p_smooth_sem"}.issubset(ps_aggregate.columns):
        wt_smooth = wt.rename(
            columns={
                "p_smooth_mean": "p_smooth_wt_mean",
                "p_smooth_sem": "p_smooth_wt_sem",
            }
        )
        merged = merged.merge(
            wt_smooth[["s", "p_smooth_wt_mean", "p_smooth_wt_sem"]],
            on="s",
            how="left",
        )
        merged["delta_p_smooth"] = merged["p_smooth_mean"] - merged["p_smooth_wt_mean"]
        merged["delta_p_smooth_sem"] = np.sqrt(
            merged["p_smooth_sem"] ** 2 + merged["p_smooth_wt_sem"] ** 2
        )
    return merged


def _aggregate_contact_maps(runs: tuple[RunAnalysis, ...]) -> dict[str, np.ndarray]:
    means: dict[str, np.ndarray] = {}
    for condition, grouped in _group_runs(runs).items():
        maps = np.stack([run.contact_map for run in grouped], axis=0)
        means[condition] = maps.mean(axis=0)
    return means


def _compare_contact_eigs(
    runs: tuple[RunAnalysis, ...],
    wt_condition: str,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    by_condition: dict[str, list[pd.DataFrame]] = {}
    for run in runs:
        path = run.analysis_dir / "contact_eigs.csv"
        if path.is_file():
            by_condition.setdefault(run.condition, []).append(pd.read_csv(path))
    if wt_condition not in by_condition:
        return None
    condition_metadata = _condition_metadata(runs)
    means = {
        condition: _mean_ev_tables(tables)
        for condition, tables in by_condition.items()
        if tables
    }
    wt = means[wt_condition]
    summary_rows: list[dict[str, float | int | str]] = []
    delta_rows: list[pd.DataFrame] = []
    for condition, table in means.items():
        if condition == wt_condition:
            continue
        correlation = ev1_correlation(wt, table)
        merged = wt[["residue_index", "EV1"]].merge(
            table[["residue_index", "EV1"]],
            on="residue_index",
            suffixes=("_wt", ""),
        )
        merged["delta_EV1"] = merged["EV1"] - merged["EV1_wt"]
        merged["condition"] = condition
        delta_rows.append(merged[["condition", "residue_index", "EV1_wt", "EV1", "delta_EV1"]])
        summary_rows.append(
            {
                "condition": condition,
                "wt_condition": wt_condition,
                **condition_metadata.get(condition, {}),
                "EV1_correlation": correlation,
                "n_residues": int(len(merged)),
            }
        )
    if not summary_rows:
        return None
    return pd.DataFrame(summary_rows), pd.concat(delta_rows, ignore_index=True)


def _mean_ev_tables(tables: list[pd.DataFrame]) -> pd.DataFrame:
    combined = pd.concat(tables, ignore_index=True)
    value_columns = [column for column in combined.columns if column.startswith("EV")]
    grouped = combined.groupby("residue_index", as_index=False)[value_columns].mean()
    names = combined[["residue_index", "residue_name"]].drop_duplicates("residue_index")
    return grouped.merge(names, on="residue_index", how="left")


def _compare_pca_centroids(runs: tuple[RunAnalysis, ...], wt_condition: str) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    condition_metadata = _condition_metadata(runs)
    for source in ("feature_pca", "contact_pca"):
        by_condition: dict[str, list[np.ndarray]] = {}
        for run in runs:
            path = run.analysis_dir / f"{source}_scores.parquet"
            if not path.is_file():
                continue
            scores = pd.read_parquet(path)
            if {"PC1", "PC2"}.issubset(scores.columns):
                by_condition.setdefault(run.condition, []).append(
                    scores[["PC1", "PC2"]].mean(axis=0).to_numpy(dtype=float)
                )
        if wt_condition not in by_condition:
            continue
        wt_centroid = np.mean(np.stack(by_condition[wt_condition], axis=0), axis=0)
        for condition, centroids in by_condition.items():
            if condition == wt_condition:
                continue
            centroid = np.mean(np.stack(centroids, axis=0), axis=0)
            delta = centroid - wt_centroid
            rows.append(
                {
                    "condition": condition,
                    "wt_condition": wt_condition,
                    "source": source,
                    **condition_metadata.get(condition, {}),
                    "delta_PC1": float(delta[0]),
                    "delta_PC2": float(delta[1]),
                    "centroid_shift": float(np.linalg.norm(delta)),
                }
            )
    return pd.DataFrame(rows)


def _condition_metadata(runs: tuple[object, ...]) -> dict[str, dict[str, float | int | str]]:
    metadata: dict[str, dict[str, float | int | str]] = {}
    for condition, grouped in _group_runs(runs).items():
        states = sorted(
            {
                str(getattr(run, "cleavage_state", ""))
                for run in grouped
                if getattr(run, "cleavage_state", "")
            }
        )
        cut_numbers: list[int] = []
        event_times: list[float] = []
        for run in grouped:
            cut_number = getattr(run, "cut_number", None)
            event_time = getattr(run, "event_time_ns", None)
            if cut_number is not None:
                cut_numbers.append(int(cut_number))
            if event_time is not None:
                event_times.append(float(event_time))
        row: dict[str, float | int | str] = {}
        if states:
            row["cleavage_state"] = ";".join(states)
        if cut_numbers:
            row["cut_number"] = int(round(float(np.mean(cut_numbers))))
        if event_times:
            row["event_time_ns"] = float(np.mean(event_times))
        metadata[condition] = row
    return metadata


def _group_runs(runs: tuple[RunAnalysis, ...]) -> dict[str, tuple[RunAnalysis, ...]]:
    groups: dict[str, list[RunAnalysis]] = {}
    for run in runs:
        groups.setdefault(run.condition, []).append(run)
    return {condition: tuple(grouped) for condition, grouped in groups.items()}


def _slug(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return clean.strip("._-") or "condition"
