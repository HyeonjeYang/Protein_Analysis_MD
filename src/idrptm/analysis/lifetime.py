"""Contact-lifetime analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from idrptm.analysis._validation import validate_lag_count


def contact_lifetime(
    contacts: ArrayLike,
    max_lag: int | None = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """Return intermittent contact autocorrelation by lag.

    ``contacts`` may have shape ``(n_frames,)``, ``(n_frames, n_pairs)``, or
    ``(n_frames, n_residues, n_residues)``. Non-zero values are treated as
    contacts. Matrix inputs are reduced to upper-triangle residue pairs.
    """

    states = _as_contact_states(contacts)
    n_frames = states.shape[0]
    max_lag_value = validate_lag_count(n_frames, max_lag)
    occupancy = states.mean()

    rows: list[dict[str, float | int]] = []
    for lag in range(max_lag_value + 1):
        if lag == 0:
            products = states * states
        else:
            products = states[:-lag] * states[lag:]
        raw = float(products.mean()) if products.size else float("nan")
        correlation = raw
        if normalize:
            correlation = raw / occupancy if occupancy > 0 else float("nan")
        rows.append(
            {
                "lag": lag,
                "correlation": correlation,
                "raw_probability": raw,
                "n_origins": int(max(n_frames - lag, 0)),
            }
        )
    return pd.DataFrame(rows, columns=["lag", "correlation", "raw_probability", "n_origins"])


def compute_contact_lifetime(
    contacts: ArrayLike,
    max_lag: int | None = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """Backward-compatible wrapper for :func:`contact_lifetime`."""

    return contact_lifetime(contacts, max_lag=max_lag, normalize=normalize)


def _as_contact_states(contacts: ArrayLike) -> np.ndarray:
    array = np.asarray(contacts)
    if array.ndim == 1:
        states = array[:, np.newaxis]
    elif array.ndim == 2:
        states = array
    elif array.ndim == 3 and array.shape[1] == array.shape[2]:
        pair_indices = np.triu_indices(array.shape[1], k=1)
        states = array[:, pair_indices[0], pair_indices[1]]
    else:
        raise ValueError(
            "contacts must have shape (n_frames,), (n_frames, n_pairs), "
            "or (n_frames, n_residues, n_residues)."
        )
    if states.shape[0] < 1:
        raise ValueError("contacts must contain at least one frame.")
    return (states != 0).astype(float)
