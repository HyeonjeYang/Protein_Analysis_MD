"""Pure analysis core for trajectory-derived observables."""

from __future__ import annotations

from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.io import TrajectoryData, load_calvados_trajectory
from idrptm.analysis.lifetime import contact_lifetime
from idrptm.analysis.msd import com_msd
from idrptm.analysis.multichain import (
    chain_com_msd,
    cluster_size_timeseries,
    com_distance_timeseries,
    inter_chain_contact_map,
    inter_protein_contact_lifetime,
    intra_chain_contact_map,
    per_chain_ree,
    per_chain_rg,
)
from idrptm.analysis.pipeline import AnalysisResult, analyze_run_directory, analyze_trajectory_data
from idrptm.analysis.ps import p_of_s
from idrptm.analysis.ree import ree_timeseries
from idrptm.analysis.rg import rg_timeseries
from idrptm.analysis.scaling import (
    FloryFit,
    fit_flory_exponent,
    internal_distance_scaling,
)

__all__ = [
    "FloryFit",
    "AnalysisResult",
    "TrajectoryData",
    "analyze_run_directory",
    "analyze_trajectory_data",
    "chain_com_msd",
    "cluster_size_timeseries",
    "com_msd",
    "com_distance_timeseries",
    "contact_lifetime",
    "contact_map_from_positions",
    "fit_flory_exponent",
    "inter_chain_contact_map",
    "inter_protein_contact_lifetime",
    "internal_distance_scaling",
    "intra_chain_contact_map",
    "load_calvados_trajectory",
    "p_of_s",
    "per_chain_ree",
    "per_chain_rg",
    "ree_timeseries",
    "rg_timeseries",
]
