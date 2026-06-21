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
    "minimal": {"observables": ["rg", "ree", "energy"]},
    "standard_idr": {
        "observables": ["rg", "ree", "contacts", "ps", "scaling", "msd", "energy"],
        "contact_cutoff_nm": 0.8,
        "min_sequence_separation": 1,
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
    },
}

REPORT_PRESETS: dict[ReportPresetName, dict[str, Any]] = {
    "minimal": {"formats": ["png"], "html": False},
    "standard": {"formats": ["png", "pdf"], "html": False},
    "publication_draft": {"formats": ["png", "pdf"], "html": True},
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
