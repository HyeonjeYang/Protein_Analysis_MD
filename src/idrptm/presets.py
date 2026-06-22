"""Preset registries for concise protein_analysis_md configs."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

SimulationPresetName = str
AnalysisPresetName = str
ReportPresetName = str


SIMULATION_PRESETS: dict[SimulationPresetName, dict[str, Any]] = {
    "smoke_single_chain": {
        "equilibration": {"total_time_ns": 0.1, "frame_interval_ns": 0.1, "save_trajectory": False},
        "production": {"total_time_ns": 0.5, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 1,
        "box_nm": [34.4, 34.4, 34.4],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
    },
    "short_single_chain": {
        "equilibration": {"total_time_ns": 1.0, "frame_interval_ns": 0.5, "save_trajectory": False},
        "production": {"total_time_ns": 5.0, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 1,
        "box_nm": [34.4, 34.4, 34.4],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
    },
    "production_single_chain": {
        "equilibration": {
            "total_time_ns": 10.0,
            "frame_interval_ns": 1.0,
            "save_trajectory": False,
        },
        "production": {"total_time_ns": 100.0, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 3,
        "box_nm": [34.4, 34.4, 34.4],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
    },
    "cleavage_smoke": {
        "equilibration": {"total_time_ns": 0.1, "frame_interval_ns": 0.1, "save_trajectory": False},
        "production": {"total_time_ns": 0.5, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 1,
        "box_nm": [34.4, 34.4, 34.4],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
    },
    "cleavage_production": {
        "equilibration": {
            "total_time_ns": 10.0,
            "frame_interval_ns": 1.0,
            "save_trajectory": False,
        },
        "production": {"total_time_ns": 100.0, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 3,
        "box_nm": [34.4, 34.4, 34.4],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
    },
    "phase_smoke": {
        "equilibration": {"total_time_ns": 0.1, "frame_interval_ns": 0.1, "save_trajectory": False},
        "production": {"total_time_ns": 0.5, "frame_interval_ns": 0.1, "save_trajectory": True},
        "replicates": 1,
        "box_nm": [50.0, 50.0, 200.0],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
        "topol": "slab",
    },
    "phase_slab_production": {
        "equilibration": {
            "total_time_ns": 20.0,
            "frame_interval_ns": 1.0,
            "save_trajectory": False,
        },
        "production": {"total_time_ns": 200.0, "frame_interval_ns": 0.2, "save_trajectory": True},
        "replicates": 3,
        "box_nm": [50.0, 50.0, 300.0],
        "platform": "CPU",
        "model": "CALVADOS2",
        "dt_ps": 0.01,
        "topol": "slab",
    },
}

ANALYSIS_PRESETS: dict[AnalysisPresetName, dict[str, Any]] = {
    "minimal": {
        "observables": ["rg", "ree", "energy"],
        "smoothing": {
            "contact_map": {"enabled": False, "visualization_only": True},
            "delta_contact_map": {"enabled": False, "visualization_only": True},
        },
        "decomposition": {"enabled": False},
        "free_energy": {"enabled": False},
    },
    "standard_idr": {
        "observables": ["rg", "ree", "contacts", "ps", "scaling", "msd", "energy"],
        "contact_cutoff_nm": 0.8,
        "min_sequence_separation": 1,
        "smoothing": {
            "ps": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "rs": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "energy": {
                "enabled": True,
                "method": "rolling",
                "window": 25,
                "visualization_only": True,
            },
            "contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
            "delta_contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
        },
        "free_energy": {"enabled": False},
        "decomposition": {
            "enabled": True,
            "coordinate_pca": {
                "enabled": False,
                "remove_com": True,
                "align": "none",
                "n_components": 5,
            },
            "contact_pca": {
                "enabled": False,
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "distance_pca": {
                "enabled": False,
                "transform": "log",
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "feature_pca": {"enabled": True, "standardize_features": True, "n_components": 5},
            "contact_eigs": {
                "enabled": True,
                "oe_method": "log_ratio",
                "eps": 1.0e-6,
                "min_sequence_separation": 2,
                "orient_by": "charge_window",
                "window": 7,
                "n_eigs": 3,
            },
            "nmf": {"enabled": False, "n_modules": 3, "random_seed": 123},
        },
    },
    "standard_cleavage": {
        "observables": [
            "rg",
            "ree",
            "contacts",
            "ps",
            "scaling",
            "msd",
            "energy",
            "fragment_metrics",
            "cleavage_events",
        ],
        "contact_cutoff_nm": 0.8,
        "min_sequence_separation": 1,
        "smoothing": {
            "ps": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "rs": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "energy": {
                "enabled": True,
                "method": "rolling",
                "window": 25,
                "visualization_only": True,
            },
            "contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
            "delta_contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
        },
        "free_energy": {"enabled": False},
        "decomposition": {
            "enabled": True,
            "coordinate_pca": {
                "enabled": False,
                "remove_com": True,
                "align": "none",
                "n_components": 5,
            },
            "contact_pca": {
                "enabled": True,
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "distance_pca": {
                "enabled": False,
                "transform": "log",
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "feature_pca": {"enabled": True, "standardize_features": True, "n_components": 5},
            "contact_eigs": {
                "enabled": True,
                "oe_method": "log_ratio",
                "eps": 1.0e-6,
                "min_sequence_separation": 2,
                "orient_by": "charge_window",
                "window": 7,
                "n_eigs": 3,
            },
            "nmf": {"enabled": False, "n_modules": 3, "random_seed": 123},
        },
    },
    "standard_phase": {
        "observables": [
            "density_profile",
            "cluster_size",
            "inter_chain_contacts",
            "com_msd",
            "energy",
        ],
        "contact_cutoff_nm": 0.8,
        "min_sequence_separation": 1,
        "smoothing": {
            "density_profile": {
                "enabled": True,
                "method": "gaussian",
                "sigma_bins": 1.0,
                "visualization_only": True,
            },
            "energy": {
                "enabled": True,
                "method": "rolling",
                "window": 25,
                "visualization_only": True,
            },
            "contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
        },
        "free_energy": {"enabled": False},
        "decomposition": {
            "enabled": True,
            "feature_pca": {"enabled": True, "standardize_features": True, "n_components": 5},
        },
    },
    "full": {
        "observables": [
            "rg",
            "ree",
            "contacts",
            "ps",
            "scaling",
            "lifetime",
            "msd",
            "energy",
            "fragment_metrics",
            "phase",
        ],
        "contact_cutoff_nm": 0.8,
        "min_sequence_separation": 1,
        "smoothing": {
            "ps": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "rs": {
                "enabled": True,
                "method": "logspace",
                "window_log10": 0.2,
                "min_points": 5,
                "robust": True,
            },
            "rg": {"enabled": True, "method": "rolling", "window": 25},
            "ree": {"enabled": True, "method": "rolling", "window": 25},
            "energy": {"enabled": True, "method": "rolling", "window": 25},
            "contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
            "delta_contact_map": {
                "enabled": False,
                "method": "gaussian",
                "sigma": 1.0,
                "visualization_only": True,
            },
        },
        "free_energy": {
            "enabled": True,
            "variables": [["Rg", "Ree"]],
            "bins": 50,
            "temperature_K": 293.0,
            "min_count": 1,
        },
        "decomposition": {
            "enabled": True,
            "coordinate_pca": {
                "enabled": True,
                "remove_com": True,
                "align": "none",
                "n_components": 5,
            },
            "contact_pca": {
                "enabled": True,
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "distance_pca": {
                "enabled": True,
                "transform": "log",
                "min_sequence_separation": 2,
                "standardize_features": True,
                "n_components": 5,
            },
            "feature_pca": {"enabled": True, "standardize_features": True, "n_components": 5},
            "contact_eigs": {
                "enabled": True,
                "oe_method": "log_ratio",
                "eps": 1.0e-6,
                "min_sequence_separation": 2,
                "orient_by": "charge_window",
                "window": 7,
                "n_eigs": 3,
            },
            "nmf": {"enabled": False, "n_modules": 3, "random_seed": 123},
        },
    },
}

REPORT_PRESETS: dict[ReportPresetName, dict[str, Any]] = {
    "minimal": {
        "formats": ["png"],
        "html": False,
        "plots": {
            "single_chain": True,
            "ptm_comparison": False,
            "cleavage": False,
            "phase": False,
            "heatmaps": True,
            "decomposition": False,
            "free_energy": False,
        },
        "smoothing": {
            "use_smoothed_ps": False,
            "use_smoothed_rs": False,
            "show_raw_points": True,
            "show_smoothed_line": False,
            "show_smoothing_metadata": False,
        },
    },
    "standard": {
        "formats": ["png"],
        "html": False,
        "plots": {
            "single_chain": True,
            "ptm_comparison": True,
            "cleavage": True,
            "phase": True,
            "heatmaps": True,
            "decomposition": True,
            "free_energy": True,
        },
        "smoothing": {
            "use_smoothed_ps": True,
            "use_smoothed_rs": True,
            "show_raw_points": True,
            "show_smoothed_line": True,
            "show_smoothing_metadata": True,
        },
    },
    "publication_draft": {
        "formats": ["png", "pdf"],
        "html": True,
        "plots": {
            "single_chain": True,
            "ptm_comparison": True,
            "cleavage": True,
            "phase": True,
            "heatmaps": True,
            "decomposition": True,
            "free_energy": True,
        },
        "smoothing": {
            "use_smoothed_ps": True,
            "use_smoothed_rs": True,
            "show_raw_points": True,
            "show_smoothed_line": True,
            "show_smoothing_metadata": True,
        },
    },
}


def simulation_preset(name: str) -> dict[str, Any]:
    """Return a copy of a simulation preset or raise a readable error."""

    return _preset(SIMULATION_PRESETS, name, "simulation")


def analysis_preset(name: str) -> dict[str, Any]:
    """Return a copy of an analysis preset or raise a readable error."""

    return _preset(ANALYSIS_PRESETS, name, "analysis")


def report_preset(name: str) -> dict[str, Any]:
    """Return a copy of a report preset or raise a readable error."""

    return _preset(REPORT_PRESETS, name, "report")


def merge_overrides(defaults: dict[str, Any], overrides: dict[str, Any] | None) -> dict[str, Any]:
    """Recursively merge user overrides onto preset defaults."""

    merged = deepcopy(defaults)
    for key, value in (overrides or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_overrides(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _preset(registry: dict[str, dict[str, Any]], name: str, kind: str) -> dict[str, Any]:
    if name not in registry:
        known = ", ".join(sorted(registry))
        raise ValueError(f"Unknown {kind} preset {name!r}. Known presets: {known}.")
    return deepcopy(registry[name])
