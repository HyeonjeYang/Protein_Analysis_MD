"""Typed configuration models for idr-ptm-md workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

PTMKind = Literal["pSer", "pThr"]
VariantMode = Literal["wt", "explicit", "single_site_scan", "all_sites"]
ExecutionBackend = Literal["local", "slurm"]
Integrator = Literal["calvados_default"]
SimulationPlatform = Literal["CPU", "CUDA"]


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


class SimulationConfig(StrictModel):
    """CALVADOS simulation timing and output controls for the MVP adapter."""

    integrator: Integrator = "calvados_default"
    dt_ps: float = 0.01
    save_every_steps: int | None = Field(7000, gt=0)
    n_frames: int | None = Field(1010, gt=0)
    total_time_ns: float | None = Field(None, gt=0)
    frame_interval_ns: float | None = Field(None, gt=0)
    total_steps: int | None = Field(None, gt=0)
    runtime_hours: float = Field(0, ge=0)
    platform: SimulationPlatform = "CPU"
    restart: str | None = "checkpoint"
    checkpoint_file: str = "restart.chk"
    random_seed: int | None = None

    @field_validator("integrator", mode="before")
    @classmethod
    def reject_custom_integrator(cls, value: object) -> object:
        """Reject integrator switching until the CALVADOS adapter supports it."""

        if value != "calvados_default":
            raise ValueError(
                "MVP supports only integrator='calvados_default'; custom dt/integrator "
                "switching requires a future CALVADOS adapter extension."
            )
        return value

    @field_validator("dt_ps", mode="before")
    @classmethod
    def reject_custom_timestep(cls, value: object) -> object:
        """Reject custom timesteps until the CALVADOS adapter supports them."""

        if float(value) != 0.01:
            raise ValueError(
                "MVP supports only dt_ps=0.01; custom dt/integrator switching "
                "requires a future CALVADOS adapter extension."
            )
        return value

    @model_validator(mode="after")
    def derive_timing(self) -> SimulationConfig:
        """Derive CALVADOS steps/frame timing from the supported input modes."""

        time_mode = self.total_time_ns is not None or self.frame_interval_ns is not None
        if time_mode:
            if self.total_time_ns is None or self.frame_interval_ns is None:
                raise ValueError(
                    "Provide both total_time_ns and frame_interval_ns for time-based "
                    "simulation timing."
                )
            save_every_steps = round(self.frame_interval_ns * 1000.0 / self.dt_ps)
            n_frames = round(self.total_time_ns / self.frame_interval_ns)
            if save_every_steps <= 0 or n_frames <= 0:
                raise ValueError("Derived save_every_steps and n_frames must be positive.")
            total_steps = save_every_steps * n_frames
            self.save_every_steps = save_every_steps
            self.n_frames = n_frames
            self.total_steps = total_steps
            self.frame_interval_ns = save_every_steps * self.dt_ps / 1000.0
            self.total_time_ns = total_steps * self.dt_ps / 1000.0
            return self

        if self.save_every_steps is None or self.n_frames is None:
            raise ValueError(
                "Provide save_every_steps and n_frames, or total_time_ns and "
                "frame_interval_ns."
            )

        total_steps = self.save_every_steps * self.n_frames
        if self.total_steps is not None and self.total_steps != total_steps:
            raise ValueError(
                "total_steps must equal save_every_steps * n_frames when all three "
                "are provided."
            )
        self.total_steps = total_steps
        self.frame_interval_ns = self.save_every_steps * self.dt_ps / 1000.0
        self.total_time_ns = total_steps * self.dt_ps / 1000.0
        return self

    def metadata(self) -> dict[str, float | int | str | None]:
        """Return raw and derived simulation timing values for metadata.json."""

        return {
            "integrator": self.integrator,
            "dt_ps": self.dt_ps,
            "save_every_steps": self.save_every_steps,
            "n_frames": self.n_frames,
            "total_time_ns": self.total_time_ns,
            "frame_interval_ns": self.frame_interval_ns,
            "total_steps": self.total_steps,
            "runtime_hours": self.runtime_hours,
            "platform": self.platform,
            "restart": self.restart,
            "checkpoint_file": self.checkpoint_file,
            "random_seed": self.random_seed,
        }


class CalvadosConfig(StrictModel):
    """External CALVADOS backend settings."""

    installation: Path | None = Field(
        None,
        description="Optional path to an existing CALVADOS checkout or environment.",
    )
    model: str = "CALVADOS2"
    residue_parameters: Path | None = None
    box_nm: list[float] = Field(default_factory=lambda: [50.0, 50.0, 50.0])
    temperature_k: float = 293.0
    ph: float = 7.4
    ionic_m: float = 0.19
    topol: str = "center"
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
    verbose: bool = True
    charge_termini: Literal["both", "N", "C", "none"] = "both"
    molecule_type: str = "protein"
    nmol: int = 1

    @model_validator(mode="after")
    def validate_box(self) -> CalvadosConfig:
        """Ensure CALVADOS box dimensions are a three-vector."""

        if len(self.box_nm) != 3:
            raise ValueError("calvados.box_nm must contain exactly three dimensions.")
        return self


class RunnerConfig(StrictModel):
    """Execution scaffold settings."""

    backend: ExecutionBackend = "local"
    work_dir: Path = Path("runs")
    dry_run: bool = True


class AnalysisConfig(StrictModel):
    """Analysis settings for trajectory observables."""

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
    contact_cutoff_nm: float = 0.8
    min_sequence_separation: int = 1
    max_lag: int | None = None
    fit_min_s: int | None = None
    fit_max_s: int | None = None


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
    if (
        config.calvados.residue_parameters is not None
        and not config.calvados.residue_parameters.is_absolute()
    ):
        calvados = config.calvados.model_copy(
            update={
                "residue_parameters": config_path.parent / config.calvados.residue_parameters
            }
        )
        config = config.model_copy(update={"calvados": calvados})
    return config
