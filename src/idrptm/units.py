"""Canonical unit definitions and analysis-output unit metadata."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

CANONICAL_UNITS: dict[str, str] = {
    "length": "nm",
    "time_step": "ps",
    "time_report": "ns",
    "energy": "kJ/mol",
    "temperature": "K",
    "ionic_strength": "M",
    "charge": "e",
    "mass": "amu",
    "rg": "nm",
    "ree": "nm",
    "internal_distance": "nm",
    "internal_distance_squared": "nm^2",
    "rs": "nm",
    "msd": "nm^2",
    "contact_probability": "dimensionless",
    "ps": "dimensionless",
}

ANALYSIS_OUTPUT_UNITS: dict[str, dict[str, str]] = {
    "timeseries_rg": {
        "frame": "index",
        "time_ps": "ps",
        "rg": CANONICAL_UNITS["rg"],
        "rg_nm_smooth": CANONICAL_UNITS["rg"],
        "rg_smoothing_method": "identifier",
    },
    "timeseries_ree": {
        "frame": "index",
        "time_ps": "ps",
        "ree": CANONICAL_UNITS["ree"],
        "ree_nm_smooth": CANONICAL_UNITS["ree"],
        "ree_smoothing_method": "identifier",
    },
    "contact_map": {
        "contact_probability": CANONICAL_UNITS["contact_probability"],
    },
    "ps": {
        "s": "residues",
        "p": CANONICAL_UNITS["ps"],
        "p_contact": CANONICAL_UNITS["ps"],
        "p_contact_smooth": CANONICAL_UNITS["ps"],
        "n_pairs": "count",
        "n_points_smooth_window": "count",
        "smoothing_method": "identifier",
        "smoothing_window_log10": "log10(residues)",
    },
    "scaling": {
        "s": "residues",
        "distance": CANONICAL_UNITS["internal_distance"],
        "mean_distance_nm": CANONICAL_UNITS["internal_distance"],
        "rms_distance_nm": CANONICAL_UNITS["internal_distance"],
        "mean_square_distance_nm2": CANONICAL_UNITS["internal_distance_squared"],
        "std_distance_nm": CANONICAL_UNITS["internal_distance"],
        "mean_distance_nm_smooth": CANONICAL_UNITS["internal_distance"],
        "n_pairs": "count",
        "n_points_smooth_window": "count",
        "smoothing_method": "identifier",
        "smoothing_window_log10": "log10(residues)",
    },
    "msd": {
        "lag": "frames",
        "msd": CANONICAL_UNITS["msd"],
        "n_origins": "count",
    },
    "contact_lifetime": {
        "lag": "frames",
        "correlation": "dimensionless",
        "raw_probability": "dimensionless",
        "n_origins": "count",
    },
    "per_chain_rg": {
        "frame": "index",
        "chain_id": "identifier",
        "rg": CANONICAL_UNITS["rg"],
    },
    "per_chain_ree": {
        "frame": "index",
        "chain_id": "identifier",
        "ree": CANONICAL_UNITS["ree"],
    },
    "intra_chain_contact_map": {
        "contact_probability": CANONICAL_UNITS["contact_probability"],
    },
    "inter_chain_contact_map": {
        "contact_probability": CANONICAL_UNITS["contact_probability"],
    },
    "com_distance": {
        "frame": "index",
        "distance": CANONICAL_UNITS["length"],
    },
    "chain_com_msd": {
        "lag": "frames",
        "msd": CANONICAL_UNITS["msd"],
        "n_origins": "count",
    },
    "cluster_size": {
        "frame": "index",
        "n_clusters": "count",
        "largest_cluster_size": "chains",
        "mean_cluster_size": "chains",
    },
    "inter_protein_contact_lifetime": {
        "lag": "frames",
        "correlation": "dimensionless",
        "raw_probability": "dimensionless",
        "n_origins": "count",
    },
    "fragment_rg": {
        "frame": "index",
        "fragment_id": "identifier",
        "rg": CANONICAL_UNITS["rg"],
    },
    "fragment_ree": {
        "frame": "index",
        "fragment_id": "identifier",
        "ree": CANONICAL_UNITS["ree"],
    },
    "inter_fragment_contact_map": {
        "contact_probability": CANONICAL_UNITS["contact_probability"],
    },
    "fragment_cluster_size": {
        "frame": "index",
        "n_clusters": "count",
        "largest_cluster_size": "fragments",
        "mean_cluster_size": "fragments",
    },
    "delta_contact_map": {
        "contact_probability_delta": "dimensionless",
    },
    "original_sequence_contact_map": {
        "contact_probability": CANONICAL_UNITS["contact_probability"],
        "residue_index": "original sequence residue",
    },
    "energy": {
        "time_ns": "ns",
        "step": "count",
        "potential_energy_kj_mol": CANONICAL_UNITS["energy"],
        "potential_energy_kj_mol_smooth": CANONICAL_UNITS["energy"],
        "kinetic_energy_kj_mol": CANONICAL_UNITS["energy"],
        "kinetic_energy_kj_mol_smooth": CANONICAL_UNITS["energy"],
        "total_energy_kj_mol": CANONICAL_UNITS["energy"],
        "total_energy_kj_mol_smooth": CANONICAL_UNITS["energy"],
        "temperature_K": CANONICAL_UNITS["temperature"],
        "temperature_K_smooth": CANONICAL_UNITS["temperature"],
        "energy_smoothing_method": "identifier",
    },
    "phase_density": {
        "z_nm": CANONICAL_UNITS["length"],
        "count": "beads",
    },
    "free_energy": {
        "x": "variable-dependent",
        "y": "variable-dependent",
        "count": "count",
        "probability": "dimensionless",
        "free_energy_kj_mol": CANONICAL_UNITS["energy"],
    },
}


def analysis_output_units(output_name: str) -> dict[str, str]:
    """Return unit metadata for a named analysis output."""

    return dict(ANALYSIS_OUTPUT_UNITS[output_name])


def units_sidecar_path(path: str | Path) -> Path:
    """Return the JSON sidecar path used for unit metadata."""

    output_path = Path(path)
    return output_path.with_name(f"{output_path.name}.units.json")


def write_units_metadata(
    path: str | Path,
    units: Mapping[str, str],
    *,
    kind: str = "analysis_output",
) -> Path:
    """Write a compact JSON unit sidecar for an analysis output."""

    sidecar = units_sidecar_path(path)
    sidecar.write_text(
        json.dumps({"kind": kind, "path": str(path), "units": dict(units)}, indent=2) + "\n",
        encoding="utf-8",
    )
    return sidecar


def summary_units(
    *,
    input_position_unit: str,
    canonical_position_unit: str = CANONICAL_UNITS["length"],
) -> dict[str, object]:
    """Build the units block written to analysis summary files."""

    return {
        "canonical": dict(CANONICAL_UNITS),
        "trajectory": {
            "input_position_unit": input_position_unit,
            "canonical_position_unit": canonical_position_unit,
        },
        "analysis_outputs": {
            name: analysis_output_units(name) for name in sorted(ANALYSIS_OUTPUT_UNITS)
        },
    }


def angstrom_to_nm(values: object) -> object:
    """Convert coordinate values from Angstrom to nanometers."""

    return values / 10.0
