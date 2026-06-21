"""Pure analysis core for trajectory-derived observables."""

from __future__ import annotations

from idrptm.analysis.cleavage import (
    fragment_cluster_size,
    fragment_resolved_ree,
    fragment_resolved_rg,
    intact_vs_cleaved_delta_contact_map,
    inter_fragment_contact_map,
    original_sequence_coordinate_contact_map,
)
from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.energy import parse_energy_log, smooth_energy_timeseries, write_energy_outputs
from idrptm.analysis.equilibration import (
    block_means,
    equilibration_diagnostics,
    running_mean,
    write_equilibration_outputs,
)
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
from idrptm.analysis.phase import density_profile_z, largest_cluster_fraction
from idrptm.analysis.pipeline import AnalysisResult, analyze_run_directory, analyze_trajectory_data
from idrptm.analysis.ps import p_of_s
from idrptm.analysis.ree import ree_timeseries
from idrptm.analysis.rg import rg_timeseries
from idrptm.analysis.scaling import (
    FloryFit,
    fit_flory_exponent,
    internal_distance_scaling,
)
from idrptm.analysis.smoothing import (
    bootstrap_smooth_ci,
    coarse_bin_curve,
    logspace_smooth_1d,
    rolling_smooth_1d,
    savgol_smooth_1d,
    smooth_contact_map,
)

__all__ = [
    "FloryFit",
    "AnalysisResult",
    "TrajectoryData",
    "analyze_run_directory",
    "analyze_trajectory_data",
    "chain_com_msd",
    "cluster_size_timeseries",
    "block_means",
    "bootstrap_smooth_ci",
    "coarse_bin_curve",
    "com_msd",
    "com_distance_timeseries",
    "contact_lifetime",
    "contact_map_from_positions",
    "density_profile_z",
    "equilibration_diagnostics",
    "fit_flory_exponent",
    "fragment_cluster_size",
    "fragment_resolved_ree",
    "fragment_resolved_rg",
    "intact_vs_cleaved_delta_contact_map",
    "inter_chain_contact_map",
    "inter_fragment_contact_map",
    "inter_protein_contact_lifetime",
    "internal_distance_scaling",
    "intra_chain_contact_map",
    "load_calvados_trajectory",
    "logspace_smooth_1d",
    "largest_cluster_fraction",
    "original_sequence_coordinate_contact_map",
    "p_of_s",
    "parse_energy_log",
    "per_chain_ree",
    "per_chain_rg",
    "ree_timeseries",
    "rolling_smooth_1d",
    "running_mean",
    "rg_timeseries",
    "savgol_smooth_1d",
    "smooth_contact_map",
    "smooth_energy_timeseries",
    "write_equilibration_outputs",
    "write_energy_outputs",
]
