"""Report-generation placeholders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReportPlan:
    """A planned report artifact."""

    output: Path
    title: str = "idr-ptm-md report"


def build_report_plan(output: str | Path, title: str = "idr-ptm-md report") -> ReportPlan:
    """Create a report plan without rendering content yet."""

    return ReportPlan(output=Path(output), title=title)
