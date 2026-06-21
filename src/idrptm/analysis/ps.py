"""P(s) analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike


def p_of_s(contact_map: ArrayLike) -> pd.DataFrame:
    """Return contact probability as a function of sequence separation.

    Parameters
    ----------
    contact_map
        Square matrix of contact probabilities or binary contacts.
    """

    contacts = np.asarray(contact_map, dtype=float)
    if contacts.ndim != 2 or contacts.shape[0] != contacts.shape[1]:
        raise ValueError("contact_map must be a square matrix.")

    rows: list[dict[str, float | int]] = []
    n_residues = contacts.shape[0]
    for separation in range(1, n_residues):
        values = np.diag(contacts, k=separation)
        rows.append(
            {
                "s": separation,
                "p": float(values.mean()) if values.size else float("nan"),
                "p_contact": float(values.mean()) if values.size else float("nan"),
                "n_pairs": int(values.size),
            }
        )
    return pd.DataFrame(rows, columns=["s", "p", "p_contact", "n_pairs"])


def compute_ps(contact_map: ArrayLike) -> pd.DataFrame:
    """Backward-compatible wrapper for :func:`p_of_s`."""

    return p_of_s(contact_map)
