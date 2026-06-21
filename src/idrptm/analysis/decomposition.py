"""Exploratory PCA and contact-environment decomposition analyses."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances

AlignMode = Literal["none", "kabsch"]
OEMethod = Literal["log_ratio", "difference"]
DistanceTransform = Literal["raw", "log"]


@dataclass(frozen=True)
class PCAResult:
    """In-memory PCA result."""

    scores: pd.DataFrame
    loadings: np.ndarray
    explained_variance: pd.DataFrame
    feature_names: tuple[str, ...]


@dataclass(frozen=True)
class ContactEigenResult:
    """In-memory contact-environment eigendecomposition result."""

    contact_eigs: pd.DataFrame
    observed_expected: np.ndarray
    correlation: np.ndarray
    eigenvalues: pd.DataFrame
    orientation_metadata: dict[str, object]


def coordinate_pca(
    positions: ArrayLike,
    *,
    time_ps: ArrayLike | None = None,
    remove_com: bool = True,
    align: AlignMode = "none",
    bead_indices: ArrayLike | None = None,
    n_components: int = 5,
) -> tuple[PCAResult, pd.DataFrame]:
    """Run Cartesian coordinate PCA / essential dynamics on a trajectory."""

    trajectory = as_position_trajectory(positions)
    selected = _select_beads(trajectory, bead_indices)
    if remove_com:
        selected = selected - selected.mean(axis=1, keepdims=True)
    if align == "kabsch":
        selected = _align_kabsch(selected)
    elif align != "none":
        raise ValueError("align must be 'none' or 'kabsch'.")

    matrix = selected.reshape(selected.shape[0], -1)
    result = _pca(matrix, n_components=n_components, standardize=False)
    scores = _score_table(result.scores, time_ps=time_ps)
    result = PCAResult(
        scores=scores,
        loadings=result.loadings,
        explained_variance=result.explained_variance,
        feature_names=result.feature_names,
    )
    return result, representative_frames(scores)


def contact_pca(
    contact_trajectory: ArrayLike,
    *,
    use_upper_triangle: bool = True,
    min_sequence_separation: int = 2,
    standardize_features: bool = True,
    n_components: int = 5,
) -> tuple[PCAResult, pd.DataFrame]:
    """Run PCA on flattened per-frame contact maps."""

    contacts = _as_map_trajectory(contact_trajectory)
    matrix, pair_table = _flatten_maps(
        contacts,
        use_upper_triangle=use_upper_triangle,
        min_sequence_separation=min_sequence_separation,
    )
    result = _pca(matrix, n_components=n_components, standardize=standardize_features)
    result = PCAResult(
        scores=_score_table(result.scores),
        loadings=result.loadings,
        explained_variance=result.explained_variance,
        feature_names=tuple(_pair_name(row) for row in pair_table.itertuples()),
    )
    return result, top_loading_contacts(result.loadings, pair_table)


def distance_pca(
    distance_trajectory: ArrayLike,
    *,
    transform: DistanceTransform = "log",
    min_sequence_separation: int = 2,
    standardize_features: bool = True,
    n_components: int = 5,
) -> tuple[PCAResult, pd.DataFrame]:
    """Run PCA on flattened per-frame residue distance maps."""

    distances = _as_map_trajectory(distance_trajectory)
    if transform == "log":
        distances = np.log(np.maximum(distances, 1.0e-12))
    elif transform != "raw":
        raise ValueError("transform must be 'raw' or 'log'.")
    matrix, pair_table = _flatten_maps(
        distances,
        use_upper_triangle=True,
        min_sequence_separation=min_sequence_separation,
    )
    result = _pca(matrix, n_components=n_components, standardize=standardize_features)
    result = PCAResult(
        scores=_score_table(result.scores),
        loadings=result.loadings,
        explained_variance=result.explained_variance,
        feature_names=tuple(_pair_name(row) for row in pair_table.itertuples()),
    )
    return result, pair_table


def feature_pca(
    features: pd.DataFrame,
    *,
    standardize_features: bool = True,
    n_components: int = 5,
) -> PCAResult:
    """Run PCA on available per-frame scalar features."""

    feature_columns = [
        column
        for column in features.columns
        if column not in {"frame", "time_ps", "time_ns"}
        and pd.api.types.is_numeric_dtype(features[column])
    ]
    if not feature_columns:
        raise ValueError("feature_pca requires at least one numeric feature column.")
    table = features[feature_columns].replace([np.inf, -np.inf], np.nan)
    table = table.fillna(table.mean(numeric_only=True)).fillna(0.0)
    result = _pca(
        table.to_numpy(dtype=float),
        n_components=n_components,
        standardize=standardize_features,
        feature_names=feature_columns,
    )
    return PCAResult(
        scores=_score_table(result.scores, time_ps=features.get("time_ps")),
        loadings=result.loadings,
        explained_variance=result.explained_variance,
        feature_names=tuple(feature_columns),
    )


def contact_eigendecomposition(
    contact_map: ArrayLike,
    *,
    sequence: str | None = None,
    residue_names: list[str] | tuple[str, ...] | None = None,
    oe_method: OEMethod = "log_ratio",
    eps: float = 1.0e-6,
    min_sequence_separation: int = 2,
    orient_by: str = "charge_window",
    window: int = 7,
    n_eigs: int = 3,
) -> ContactEigenResult:
    """Compute contact-environment eigenvectors from an average contact map."""

    matrix = np.asarray(contact_map, dtype=float)
    if matrix.ndim == 3:
        matrix = matrix.mean(axis=0)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("contact_map must be square or a trajectory of square maps.")
    if eps <= 0:
        raise ValueError("eps must be positive.")
    if n_eigs <= 0:
        raise ValueError("n_eigs must be positive.")
    expected = _expected_by_separation(matrix)
    observed_expected = _observed_expected(
        matrix,
        expected,
        oe_method=oe_method,
        eps=eps,
        min_sequence_separation=min_sequence_separation,
    )
    correlation = _row_correlation(observed_expected)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    eigenvectors = eigenvectors[:, : min(n_eigs, eigenvectors.shape[1])]
    eigenvalues = eigenvalues[: eigenvectors.shape[1]]
    eigenvectors, orientation = _orient_eigenvectors(
        eigenvectors,
        sequence=sequence,
        residue_names=residue_names,
        orient_by=orient_by,
        window=window,
    )
    residue_labels = _residue_labels(matrix.shape[0], sequence, residue_names)
    eigs = pd.DataFrame(
        {
            "residue_index": np.arange(1, matrix.shape[0] + 1, dtype=int),
            "residue_name": residue_labels,
        }
    )
    for index in range(eigenvectors.shape[1]):
        eigs[f"EV{index + 1}"] = eigenvectors[:, index]
    eigen_table = pd.DataFrame(
        {
            "component": [f"EV{index + 1}" for index in range(len(eigenvalues))],
            "eigenvalue": eigenvalues,
            "explained_fraction": eigenvalues / np.sum(np.abs(eigenvalues))
            if np.sum(np.abs(eigenvalues)) > 0
            else np.zeros_like(eigenvalues),
        }
    )
    metadata = {
        "analysis_label": "contact-environment eigenvectors",
        "not_chromosome_compartments": True,
        "oe_method": oe_method,
        "eps": eps,
        "min_sequence_separation": min_sequence_separation,
        "orient_by": orient_by,
        "window": window,
        "orientation": orientation,
    }
    return ContactEigenResult(
        contact_eigs=eigs,
        observed_expected=observed_expected,
        correlation=correlation,
        eigenvalues=eigen_table,
        orientation_metadata=metadata,
    )


def nmf_contact_modules(
    contact_map: ArrayLike,
    *,
    n_modules: int = 3,
    min_sequence_separation: int = 2,
    random_seed: int = 123,
    n_iter: int = 300,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """Discover non-negative contact modules with a small multiplicative-update NMF."""

    matrix = np.asarray(contact_map, dtype=float)
    if matrix.ndim == 3:
        matrix = matrix.mean(axis=0)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("contact_map must be square or a trajectory of square maps.")
    if n_modules <= 0:
        raise ValueError("n_modules must be positive.")
    data = np.maximum(matrix.copy(), 0.0)
    pair_i, pair_j = np.indices(data.shape)
    data[np.abs(pair_i - pair_j) < min_sequence_separation] = 0.0
    rng = np.random.default_rng(random_seed)
    w = rng.random((data.shape[0], n_modules)) + 1.0e-6
    h = rng.random((n_modules, data.shape[1])) + 1.0e-6
    errors: list[dict[str, float | int]] = []
    for iteration in range(n_iter):
        reconstruction = w @ h + 1.0e-12
        h *= (w.T @ (data / reconstruction)) / np.maximum(w.sum(axis=0)[:, None], 1.0e-12)
        reconstruction = w @ h + 1.0e-12
        w *= ((data / reconstruction) @ h.T) / np.maximum(h.sum(axis=1)[None, :], 1.0e-12)
        if iteration in {0, n_iter - 1} or (iteration + 1) % 50 == 0:
            error = float(np.linalg.norm(data - w @ h))
            errors.append({"iteration": iteration + 1, "reconstruction_error": error})
    weights = pd.DataFrame({"residue_index": np.arange(1, data.shape[0] + 1)})
    for module in range(n_modules):
        weights[f"module_{module + 1}"] = w[:, module]
    pair_weights = np.stack([np.outer(w[:, index], h[index, :]) for index in range(n_modules)])
    return weights, pair_weights, pd.DataFrame(errors)


def contact_map_trajectory(
    positions: ArrayLike,
    *,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> np.ndarray:
    """Build a binary contact-map trajectory from coordinates."""

    trajectory = as_position_trajectory(positions)
    maps = np.zeros((trajectory.shape[0], trajectory.shape[1], trajectory.shape[1]), dtype=float)
    for frame_index, frame in enumerate(trajectory):
        distances = pairwise_distances(frame)
        contacts = (distances <= cutoff).astype(float)
        pair_i, pair_j = np.indices(contacts.shape)
        contacts[np.abs(pair_i - pair_j) < min_sequence_separation] = 0.0
        np.fill_diagonal(contacts, 0.0)
        maps[frame_index] = contacts
    return maps


def distance_map_trajectory(positions: ArrayLike) -> np.ndarray:
    """Build a residue distance-map trajectory from coordinates."""

    trajectory = as_position_trajectory(positions)
    return np.stack([pairwise_distances(frame) for frame in trajectory], axis=0)


def total_contacts_features(
    contact_maps: np.ndarray,
    *,
    min_long_range_separation: int = 12,
) -> pd.DataFrame:
    """Return simple contact-count features per frame."""

    n_residues = contact_maps.shape[1]
    upper_i, upper_j = np.triu_indices(n_residues, k=1)
    long_i, long_j = np.triu_indices(n_residues, k=min_long_range_separation)
    return pd.DataFrame(
        {
            "total_contacts": contact_maps[:, upper_i, upper_j].sum(axis=1),
            "long_range_contacts": contact_maps[:, long_i, long_j].sum(axis=1)
            if len(long_i)
            else np.zeros(contact_maps.shape[0]),
        }
    )


def representative_frames(scores: pd.DataFrame) -> pd.DataFrame:
    """Pick min/max/near-zero representative frames for each PC."""

    rows: list[dict[str, float | int | str]] = []
    for column in [name for name in scores.columns if name.startswith("PC")]:
        values = scores[column].to_numpy(dtype=float)
        for label, index in (
            ("minimum", int(np.nanargmin(values))),
            ("maximum", int(np.nanargmax(values))),
            ("near_zero", int(np.nanargmin(np.abs(values)))),
        ):
            rows.append(
                {
                    "component": column,
                    "representative": label,
                    "frame": int(scores.loc[index, "frame"]),
                    "score": float(values[index]),
                }
            )
    return pd.DataFrame(rows)


def top_loading_contacts(
    loadings: np.ndarray,
    pair_table: pd.DataFrame,
    *,
    n_top: int = 20,
) -> pd.DataFrame:
    """Return top positive and negative contacts for PC loadings."""

    rows: list[pd.DataFrame] = []
    for component_index, loading in enumerate(loadings, start=1):
        count = min(n_top, loading.size)
        positive = np.argsort(loading)[-count:][::-1]
        negative = np.argsort(loading)[:count]
        for label, indices in (("positive", positive), ("negative", negative)):
            table = pair_table.iloc[indices].copy()
            table["component"] = f"PC{component_index}"
            table["sign"] = label
            table["loading"] = loading[indices]
            rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def ev1_correlation(reference: pd.DataFrame, variant: pd.DataFrame) -> float:
    """Compute EV1 correlation between two contact-eigen tables."""

    merged = reference[["residue_index", "EV1"]].merge(
        variant[["residue_index", "EV1"]],
        on="residue_index",
        suffixes=("_reference", "_variant"),
    )
    if len(merged) < 2:
        return float("nan")
    return float(np.corrcoef(merged["EV1_reference"], merged["EV1_variant"])[0, 1])


def run_decomposition_analysis(
    *,
    output_dir: str | Path,
    positions: ArrayLike,
    frame_table: pd.DataFrame,
    contact_map: ArrayLike,
    analysis_config: object,
    rg: ArrayLike | None = None,
    ree: ArrayLike | None = None,
    sequence: str | None = None,
) -> dict[str, Path]:
    """Run enabled decomposition analyses and write output files."""

    config = _enabled_decomposition_config(analysis_config)
    if not config:
        return {}
    root = Path(output_dir)
    outputs: dict[str, Path] = {}
    positions_array = as_position_trajectory(positions)
    contact_maps: np.ndarray | None = None
    distance_maps: np.ndarray | None = None

    coordinate_config = _section(config, "coordinate_pca")
    if coordinate_config:
        result, frames = coordinate_pca(
            positions_array,
            time_ps=frame_table.get("time_ps"),
            remove_com=bool(coordinate_config.get("remove_com", True)),
            align=str(coordinate_config.get("align", "none")),
            n_components=int(coordinate_config.get("n_components", 5)),
        )
        outputs.update(_write_pca_result(root, "pca", result))
        variance_alias = root / "explained_variance.csv"
        result.explained_variance.to_csv(variance_alias, index=False)
        outputs["explained_variance"] = variance_alias
        path = root / "representative_frames.csv"
        frames.to_csv(path, index=False)
        outputs["representative_frames"] = path

    contact_config = _section(config, "contact_pca")
    if contact_config:
        if contact_maps is None:
            contact_maps = contact_map_trajectory(
                positions_array,
                cutoff=float(getattr(analysis_config, "contact_cutoff_nm", 0.8)),
                min_sequence_separation=int(contact_config.get("min_sequence_separation", 2)),
            )
        result, top_contacts = contact_pca(
            contact_maps,
            min_sequence_separation=int(contact_config.get("min_sequence_separation", 2)),
            standardize_features=bool(contact_config.get("standardize_features", True)),
            n_components=int(contact_config.get("n_components", 5)),
        )
        outputs.update(_write_pca_result(root, "contact_pca", result))
        path = root / "top_loading_contacts.csv"
        top_contacts.to_csv(path, index=False)
        outputs["top_loading_contacts"] = path

    distance_config = _section(config, "distance_pca")
    if distance_config:
        distance_maps = (
            distance_maps
            if distance_maps is not None
            else distance_map_trajectory(positions_array)
        )
        result, _ = distance_pca(
            distance_maps,
            transform=str(distance_config.get("transform", "log")),
            min_sequence_separation=int(distance_config.get("min_sequence_separation", 2)),
            standardize_features=bool(distance_config.get("standardize_features", True)),
            n_components=int(distance_config.get("n_components", 5)),
        )
        outputs.update(_write_pca_result(root, "distance_pca", result))

    feature_config = _section(config, "feature_pca")
    if feature_config:
        contact_maps = contact_maps if contact_maps is not None else contact_map_trajectory(
            positions_array,
            cutoff=float(getattr(analysis_config, "contact_cutoff_nm", 0.8)),
            min_sequence_separation=1,
        )
        features = frame_table.copy()
        if rg is not None:
            features["rg"] = np.asarray(rg, dtype=float)
        if ree is not None:
            features["ree"] = np.asarray(ree, dtype=float)
        features = pd.concat([features, total_contacts_features(contact_maps)], axis=1)
        result = feature_pca(
            features,
            standardize_features=bool(feature_config.get("standardize_features", True)),
            n_components=int(feature_config.get("n_components", 5)),
        )
        outputs.update(_write_pca_result(root, "feature_pca", result, loadings_csv=True))

    eig_config = _section(config, "contact_eigs")
    if eig_config:
        result = contact_eigendecomposition(
            contact_map,
            sequence=sequence,
            oe_method=str(eig_config.get("oe_method", "log_ratio")),
            eps=float(eig_config.get("eps", 1.0e-6)),
            min_sequence_separation=int(eig_config.get("min_sequence_separation", 2)),
            orient_by=str(eig_config.get("orient_by", "charge_window")),
            window=int(eig_config.get("window", 7)),
            n_eigs=int(eig_config.get("n_eigs", 3)),
        )
        outputs.update(_write_contact_eigs(root, result))

    nmf_config = _section(config, "nmf")
    if nmf_config:
        weights, pair_weights, errors = nmf_contact_modules(
            contact_map,
            n_modules=int(nmf_config.get("n_modules", 3)),
            min_sequence_separation=int(nmf_config.get("min_sequence_separation", 2)),
            random_seed=int(nmf_config.get("random_seed", 123)),
        )
        weights_path = root / "nmf_residue_weights.csv"
        pair_path = root / "nmf_pair_weights.npy"
        error_path = root / "nmf_reconstruction_error.csv"
        weights.to_csv(weights_path, index=False)
        np.save(pair_path, pair_weights)
        errors.to_csv(error_path, index=False)
        outputs |= {
            "nmf_residue_weights": weights_path,
            "nmf_pair_weights": pair_path,
            "nmf_reconstruction_error": error_path,
        }

    metadata_path = root / "decomposition_metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "exploratory": True,
                "label": "contact-environment eigenvectors and PCA analyses",
                "not_chromosome_compartments": True,
                "config": config,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["decomposition_metadata"] = metadata_path
    return outputs


def _pca(
    matrix: np.ndarray,
    *,
    n_components: int,
    standardize: bool,
    feature_names: list[str] | tuple[str, ...] | None = None,
) -> PCAResult:
    data = np.asarray(matrix, dtype=float)
    if data.ndim != 2:
        raise ValueError("PCA input must be a two-dimensional matrix.")
    if data.shape[0] < 1 or data.shape[1] < 1:
        raise ValueError("PCA input must contain at least one sample and one feature.")
    if n_components <= 0:
        raise ValueError("n_components must be positive.")
    center = np.nanmean(data, axis=0)
    centered = np.nan_to_num(data - center, nan=0.0)
    if standardize:
        scale = np.nanstd(centered, axis=0)
        scale[scale == 0] = 1.0
        centered = centered / scale
    component_count = min(n_components, data.shape[0], data.shape[1])
    u, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :component_count] * singular_values[:component_count]
    loadings = vt[:component_count]
    denominator = max(data.shape[0] - 1, 1)
    variance = singular_values[:component_count] ** 2 / denominator
    total_variance = np.sum(singular_values**2 / denominator)
    explained = pd.DataFrame(
        {
            "component": [f"PC{index + 1}" for index in range(component_count)],
            "explained_variance": variance,
            "explained_variance_ratio": variance / total_variance
            if total_variance > 0
            else np.zeros_like(variance),
        }
    )
    score_table = pd.DataFrame(
        {f"PC{index + 1}": scores[:, index] for index in range(component_count)}
    )
    names = tuple(feature_names or [f"feature_{index}" for index in range(data.shape[1])])
    return PCAResult(
        scores=score_table,
        loadings=loadings,
        explained_variance=explained,
        feature_names=names,
    )


def _select_beads(positions: np.ndarray, bead_indices: ArrayLike | None) -> np.ndarray:
    if bead_indices is None:
        return positions.copy()
    indices = np.asarray(bead_indices, dtype=int)
    return positions[:, indices, :].copy()


def _align_kabsch(positions: np.ndarray) -> np.ndarray:
    reference = positions[0]
    aligned = np.empty_like(positions)
    aligned[0] = reference
    for index, frame in enumerate(positions[1:], start=1):
        covariance = frame.T @ reference
        u, _, vt = np.linalg.svd(covariance)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        aligned[index] = frame @ rotation
    return aligned


def _score_table(scores: pd.DataFrame, *, time_ps: ArrayLike | None = None) -> pd.DataFrame:
    table = scores.copy()
    table.insert(0, "frame", np.arange(len(table), dtype=int))
    if time_ps is not None:
        values = np.asarray(time_ps, dtype=float)
        if values.shape[0] == len(table):
            table.insert(1, "time_ns", values / 1000.0)
    return table


def _as_map_trajectory(values: ArrayLike) -> np.ndarray:
    maps = np.asarray(values, dtype=float)
    if maps.ndim != 3 or maps.shape[1] != maps.shape[2]:
        raise ValueError("map trajectory must have shape (frames, residues, residues).")
    return maps


def _flatten_maps(
    maps: np.ndarray,
    *,
    use_upper_triangle: bool,
    min_sequence_separation: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    n_residues = maps.shape[1]
    if use_upper_triangle:
        pair_i, pair_j = np.triu_indices(n_residues, k=min_sequence_separation)
    else:
        pair_i, pair_j = np.where(
            np.abs(np.subtract.outer(np.arange(n_residues), np.arange(n_residues)))
            >= min_sequence_separation
        )
    pair_table = pd.DataFrame(
        {
            "residue_i": pair_i + 1,
            "residue_j": pair_j + 1,
            "sequence_separation": np.abs(pair_j - pair_i),
        }
    )
    return maps[:, pair_i, pair_j], pair_table


def _pair_name(row: object) -> str:
    return f"{int(row.residue_i)}-{int(row.residue_j)}"


def _expected_by_separation(matrix: np.ndarray) -> np.ndarray:
    expected = np.zeros_like(matrix, dtype=float)
    for separation in range(matrix.shape[0]):
        if separation == 0:
            value = 0.0
        else:
            values = np.diag(matrix, k=separation)
            value = float(np.nanmean(values)) if values.size else 0.0
        residue_indices = np.arange(matrix.shape[0])
        separations = np.abs(np.subtract.outer(residue_indices, residue_indices))
        i, j = np.where(separations == separation)
        expected[i, j] = value
    return expected


def _observed_expected(
    matrix: np.ndarray,
    expected: np.ndarray,
    *,
    oe_method: OEMethod,
    eps: float,
    min_sequence_separation: int,
) -> np.ndarray:
    if oe_method == "log_ratio":
        oe = np.log((matrix + eps) / (expected + eps))
    elif oe_method == "difference":
        oe = matrix - expected
    else:
        raise ValueError("oe_method must be 'log_ratio' or 'difference'.")
    pair_i, pair_j = np.indices(matrix.shape)
    oe[np.abs(pair_i - pair_j) < min_sequence_separation] = np.nan
    return oe


def _row_correlation(observed_expected: np.ndarray) -> np.ndarray:
    values = observed_expected.copy()
    row_means = np.nanmean(values, axis=1)
    row_means = np.nan_to_num(row_means, nan=0.0)
    missing = ~np.isfinite(values)
    values[missing] = np.take(row_means, np.where(missing)[0])
    values = values - values.mean(axis=1, keepdims=True)
    row_std = values.std(axis=1, keepdims=True)
    row_std[row_std == 0] = 1.0
    values = values / row_std
    correlation = (values @ values.T) / max(values.shape[1], 1)
    return np.nan_to_num(correlation, nan=0.0)


def _orient_eigenvectors(
    eigenvectors: np.ndarray,
    *,
    sequence: str | None,
    residue_names: list[str] | tuple[str, ...] | None,
    orient_by: str,
    window: int,
) -> tuple[np.ndarray, dict[str, object]]:
    feature = _orientation_feature(
        eigenvectors.shape[0],
        sequence=sequence,
        residue_names=residue_names,
        orient_by=orient_by,
        window=window,
    )
    oriented = eigenvectors.copy()
    flips: dict[str, bool] = {}
    for index in range(oriented.shape[1]):
        vector = oriented[:, index]
        use_feature = feature is not None and np.std(feature) > 0 and np.std(vector) > 0
        if use_feature:
            score = float(np.corrcoef(vector, feature)[0, 1])
            flip = bool(np.isfinite(score) and score < 0)
        else:
            score = float(np.sum(vector[: max(1, len(vector) // 2)]))
            flip = score < 0
        if flip:
            oriented[:, index] = -vector
        flips[f"EV{index + 1}"] = flip
    return oriented, {"feature_available": feature is not None, "flipped": flips}


def _orientation_feature(
    length: int,
    *,
    sequence: str | None,
    residue_names: list[str] | tuple[str, ...] | None,
    orient_by: str,
    window: int,
) -> np.ndarray | None:
    residues = sequence or "".join(residue_names or [])
    if len(residues) != length:
        return None
    if orient_by == "charge_window":
        values = np.array([_residue_charge(residue) for residue in residues], dtype=float)
    elif orient_by == "aromaticity":
        values = np.array([1.0 if residue in {"F", "Y", "W"} else 0.0 for residue in residues])
    elif orient_by == "hydropathy":
        values = np.array([_hydropathy(residue) for residue in residues], dtype=float)
    else:
        return None
    return _rolling_mean(values, window=max(1, window))


def _residue_charge(residue: str) -> float:
    charges = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1, "B": -2.0, "O": -2.0}
    return charges.get(residue.upper(), 0.0)


def _hydropathy(residue: str) -> float:
    values = {
        "I": 4.5,
        "V": 4.2,
        "L": 3.8,
        "F": 2.8,
        "C": 2.5,
        "M": 1.9,
        "A": 1.8,
        "G": -0.4,
        "T": -0.7,
        "S": -0.8,
        "W": -0.9,
        "Y": -1.3,
        "P": -1.6,
        "H": -3.2,
        "E": -3.5,
        "Q": -3.5,
        "D": -3.5,
        "N": -3.5,
        "K": -3.9,
        "R": -4.5,
    }
    return values.get(residue.upper(), 0.0)


def _rolling_mean(values: np.ndarray, *, window: int) -> np.ndarray:
    window = min(window, len(values))
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")


def _residue_labels(
    length: int,
    sequence: str | None,
    residue_names: list[str] | tuple[str, ...] | None,
) -> list[str]:
    if sequence and len(sequence) == length:
        return list(sequence)
    if residue_names and len(residue_names) == length:
        return list(residue_names)
    return ["X"] * length


def _enabled_decomposition_config(analysis_config: object) -> dict[str, object]:
    config = getattr(analysis_config, "decomposition", {}) or {}
    if not isinstance(config, dict) or not config.get("enabled", False):
        return {}
    return dict(config)


def _section(config: dict[str, object], key: str) -> dict[str, object]:
    value = config.get(key, {})
    if not isinstance(value, dict) or not value.get("enabled", False):
        return {}
    return dict(value)


def _write_pca_result(
    root: Path,
    prefix: str,
    result: PCAResult,
    *,
    loadings_csv: bool = False,
) -> dict[str, Path]:
    score_path = root / f"{prefix}_scores.parquet"
    loading_path = root / f"{prefix}_loadings.{'csv' if loadings_csv else 'npy'}"
    variance_path = root / f"{prefix}_explained_variance.csv"
    result.scores.to_parquet(score_path, index=False)
    if loadings_csv:
        table = pd.DataFrame(
            result.loadings.T,
            columns=[f"PC{index + 1}" for index in range(result.loadings.shape[0])],
        )
        table.insert(0, "feature", result.feature_names)
        table.to_csv(loading_path, index=False)
    else:
        np.save(loading_path, result.loadings)
    result.explained_variance.to_csv(variance_path, index=False)
    return {
        f"{prefix}_scores": score_path,
        f"{prefix}_loadings": loading_path,
        f"{prefix}_explained_variance": variance_path,
    }


def _write_contact_eigs(root: Path, result: ContactEigenResult) -> dict[str, Path]:
    eig_path = root / "contact_eigs.csv"
    oe_path = root / "contact_oe.npy"
    corr_path = root / "contact_correlation.npy"
    eigen_path = root / "eigenvalues.csv"
    metadata_path = root / "ev_orientation_metadata.json"
    result.contact_eigs.to_csv(eig_path, index=False)
    np.save(oe_path, result.observed_expected)
    np.save(corr_path, result.correlation)
    result.eigenvalues.to_csv(eigen_path, index=False)
    metadata_path.write_text(
        json.dumps(result.orientation_metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "contact_eigs": eig_path,
        "contact_oe": oe_path,
        "contact_correlation": corr_path,
        "eigenvalues": eigen_path,
        "ev_orientation_metadata": metadata_path,
    }
