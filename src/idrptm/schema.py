"""Typed configuration models for idr-ptm-md workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

PTMKind = Literal["pSer", "pThr"]
VariantMode = Literal["wt", "explicit", "single_site_scan", "all_sites"]
ExecutionBackend = Literal["local", "slurm"]


class StrictModel(BaseModel):
    """Base model that rejects unknown configuration keys."""

    model_config = ConfigDict(extra="forbid")


class SequenceConfig(StrictModel):
    """Sequence input for a single-chain MVP workflow."""

    name: str = Field(..., description="Human-readable sequence identifier.")
    sequence: str | None = Field(None, description="Inline one-letter amino-acid sequence.")
    fasta: Path | None = Field(None, description="Optional FASTA file path.")

    @model_validator(mode="after")
    def require_exactly_one_sequence_source(self) -> SequenceConfig:
        """Require either an inline sequence or a FASTA file, not both."""

        if bool(self.sequence) == bool(self.fasta):
            raise ValueError("Provide exactly one sequence source: 'sequence' or 'fasta'.")
        return self


class PTMSite(StrictModel):
    """One phosphorylatable residue position using one-based indexing."""

    position: int = Field(..., ge=1)
    residue: Literal["S", "T"]
    ptm: PTMKind

    @property
    def zero_based_index(self) -> int:
        """Internal zero-based index corresponding to ``position``."""

        return self.position - 1

    @model_validator(mode="after")
    def require_ptm_residue_match(self) -> PTMSite:
        """Ensure the configured PTM is chemically compatible with the residue."""

        expected_residue = {"pSer": "S", "pThr": "T"}[self.ptm]
        if self.residue != expected_residue:
            raise ValueError(
                f"{self.ptm} can only be configured on {expected_residue}, "
                f"not {self.residue}."
            )
        return self


class PTMConfig(StrictModel):
    """PTM design settings for the MVP phosphorylation scope."""

    mode: VariantMode = "wt"
    sites: list[PTMSite] = Field(default_factory=list)
    include_wt: bool = True

    @model_validator(mode="after")
    def validate_site_requirements(self) -> PTMConfig:
        """Validate mode-specific PTM site requirements."""

        if self.mode == "explicit" and not self.sites:
            raise ValueError("PTM mode 'explicit' requires at least one site.")
        positions = [site.position for site in self.sites]
        if len(positions) != len(set(positions)):
            raise ValueError("PTM site positions must be unique.")
        return self


class CalvadosConfig(StrictModel):
    """External CALVADOS backend settings."""

    installation: Path | None = Field(
        None,
        description="Optional path to an existing CALVADOS checkout or environment.",
    )
    model: str = "CALVADOS2"
    residue_parameters: Path | None = None
    temperature_k: float = 293.0
    ph: float = 7.4


class RunnerConfig(StrictModel):
    """Execution scaffold settings."""

    backend: ExecutionBackend = "local"
    work_dir: Path = Path("runs")
    dry_run: bool = True


class AnalysisConfig(StrictModel):
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


class WorkflowConfig(StrictModel):
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
    config = WorkflowConfig.model_validate(data)
    if config.sequence.fasta is not None and not config.sequence.fasta.is_absolute():
        sequence = config.sequence.model_copy(
            update={"fasta": config_path.parent / config.sequence.fasta}
        )
        config = config.model_copy(update={"sequence": sequence})
    return config
