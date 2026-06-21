"""Console entry point for the ``pamd`` command."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))

from idrptm.cli import app, main

__all__ = ["app", "main"]
