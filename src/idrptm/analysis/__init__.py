"""Pure analysis core for trajectory-derived observables."""

from __future__ import annotations

from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.lifetime import contact_lifetime
from idrptm.analysis.msd import com_msd
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
    "com_msd",
    "contact_lifetime",
    "contact_map_from_positions",
    "fit_flory_exponent",
    "internal_distance_scaling",
    "p_of_s",
    "ree_timeseries",
    "rg_timeseries",
]
