"""Typed configuration models for idr-ptm-md workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field

PTMKind = Literal["pSer", "pThr"]
VariantMode = Literal["wt", "explicit", "single_site_scan", "all_sites"]
ExecutionBackend = Literal["local", "slurm"]


class SequenceConfig(BaseModel):
    """Sequence input for a single-chain MVP workflow."""

    name: str = Field(..., description="Human-readable sequence identifier.")
    sequence: str | None = Field(None, description="Inline one-letter amino-acid sequence.")
    fasta: Path | None = Field(None, description="Optional FASTA file path.")


class PTMSite(BaseModel):
    """One phosphorylatable residue position using one-based indexing."""

    position: int = Field(..., ge=1)
    residue: Literal["S", "T"]
    ptm: PTMKind


class PTMConfig(BaseModel):
    """PTM design settings for the MVP phosphorylation scope."""

    mode: VariantMode = "wt"
    sites: list[PTMSite] = Field(default_factory=list)
    include_wt: bool = True


class CalvadosConfig(BaseModel):
    """External CALVADOS backend settings."""

    installation: Path | None = Field(
        None,
        description="Optional path to an existing CALVADOS checkout or environment.",
    )
    model: str = "CALVADOS2"
    residue_parameters: Path | None = None
    temperature_k: float = 293.0
    ph: float = 7.4


class RunnerConfig(BaseModel):
    """Execution scaffold settings."""

    backend: ExecutionBackend = "local"
    work_dir: Path = Path("runs")
    dry_run: bool = True


class AnalysisConfig(BaseModel):
    """Analysis outputs planned for Stage 2+."""

    observables: list[str] = Field(
        default_factory=lambda: [
            "rg",
            "ree",
            "contacts",
            "ps",
            "scaling",
            "lifetime",
            "msd",
        ]
    )


class WorkflowConfig(BaseModel):
    """Top-level idr-ptm-md workflow configuration."""

    project: str
    sequence: SequenceConfig
    ptm: PTMConfig = Field(default_factory=PTMConfig)
    calvados: CalvadosConfig = Field(default_factory=CalvadosConfig)
    runner: RunnerConfig = Field(default_factory=RunnerConfig)
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)


def load_config(path: str | Path) -> WorkflowConfig:
    """Load a workflow YAML file into a typed configuration model."""

    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if hasattr(WorkflowConfig, "model_validate"):
        return WorkflowConfig.model_validate(data)
    return WorkflowConfig.parse_obj(data)
