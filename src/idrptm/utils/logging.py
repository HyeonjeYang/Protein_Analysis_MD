"""Logging setup helpers."""

from __future__ import annotations

import logging


def get_logger(name: str = "protein_analysis_md") -> logging.Logger:
    """Return a configured package logger."""

    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger
