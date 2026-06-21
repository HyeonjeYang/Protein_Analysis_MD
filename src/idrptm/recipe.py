"""Python recipe API for defining protein_analysis_md experiments."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from idrptm.configuration import LockedConfig, compile_config, normalize_config, resolve_presets


@dataclass
class Protein:
    """Protein input definition for recipe-based configs."""

    source: str
    name: str | None = None
    sequence: str | None = None
    fasta: str | None = None
    query: str | None = None
    accession: str | None = None
    reviewed_only: bool = True
    organism: str | None = None
    interactive_select: bool = False
    region: dict[str, Any] | None = None
    charge_termini: str = "both"

    @classmethod
    def from_sequence(cls, name: str, sequence: str, **kwargs: Any) -> Protein:
        """Create a direct sequence protein."""

        return cls(source="direct", name=name, sequence=sequence, **kwargs)

    @classmethod
    def from_fasta(cls, name: str, fasta: str | Path, **kwargs: Any) -> Protein:
        """Create a FASTA-backed protein."""

        return cls(source="fasta", name=name, fasta=str(fasta), **kwargs)

    @classmethod
    def from_uniprot(
        cls,
        query: str,
        *,
        accession: str | None = None,
        reviewed_only: bool = True,
        organism: str | None = None,
        interactive_select: bool = False,
        region: dict[str, Any] | None = None,
        name: str | None = None,
        **kwargs: Any,
    ) -> Protein:
        """Create a UniProt/Swiss-Prot-backed protein."""

        return cls(
            source="uniprot",
            name=name,
            query=query,
            accession=accession,
            reviewed_only=reviewed_only,
            organism=organism,
            interactive_select=interactive_select,
            region=region,
            **kwargs,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a YAML-serializable protein block."""

        return {key: value for key, value in self.__dict__.items() if value is not None}


@dataclass
class Experiment:
    """User-facing Python experiment recipe."""

    name: str
    outdir: str | Path
    protein: Protein | None = None
    ptm_states: list[dict[str, Any]] = field(default_factory=list)
    ptm_file: str | None = None
    mutation_states: list[dict[str, Any]] = field(default_factory=list)
    cleavage: dict[str, Any] = field(default_factory=lambda: {"mode": "none"})
    protocol: dict[str, Any] = field(default_factory=lambda: {"preset": "smoke_single_chain"})
    analysis: dict[str, Any] = field(default_factory=lambda: {"preset": "standard_idr"})
    report: dict[str, Any] = field(default_factory=lambda: {"preset": "standard"})
    sweep: dict[str, Any] = field(default_factory=dict)

    def add_protein(self, protein: Protein) -> Experiment:
        """Attach the primary protein to the experiment."""

        self.protein = protein
        return self

    def use_preset(
        self,
        *,
        simulation: str | None = None,
        analysis: str | None = None,
        report: str | None = None,
    ) -> Experiment:
        """Select simulation, analysis, and report presets."""

        if simulation is not None:
            self.protocol["preset"] = simulation
        if analysis is not None:
            self.analysis["preset"] = analysis
        if report is not None:
            self.report["preset"] = report
        return self

    def add_ptm_state(
        self,
        name: str,
        modifications: list[dict[str, Any]] | None = None,
    ) -> Experiment:
        """Add a named PTM state."""

        self.ptm_states.append({"name": name, "modifications": modifications or []})
        return self

    def add_ptm_states_from_file(self, path: str | Path) -> Experiment:
        """Use PTM states from a whitespace/TSV/CSV file."""

        self.ptm_file = str(path)
        return self

    def add_mutation_state(self, name: str, mutations: list[dict[str, Any]]) -> Experiment:
        """Record a mutation state for future mutation workflows."""

        self.mutation_states.append({"name": name, "mutations": mutations})
        return self

    def add_cleavage_state(self, name: str = "intact", **kwargs: Any) -> Experiment:
        """Add one cleavage state definition."""

        self.cleavage = {"name": name, **kwargs}
        self.cleavage.setdefault("mode", "none" if name == "intact" else "manual")
        return self

    def add_environment(self, **kwargs: Any) -> Experiment:
        """Add environment overrides such as pH or ionic strength."""

        self.protocol.setdefault("environment", {}).update(kwargs)
        return self

    def add_sweep(self, **kwargs: Any) -> Experiment:
        """Record parameter sweep dimensions."""

        self.sweep.update(kwargs)
        return self

    def to_config(self) -> dict[str, Any]:
        """Return the concise user config represented by this recipe."""

        if self.protein is None:
            raise ValueError("Experiment requires a protein before it can be compiled.")
        ptm: dict[str, Any]
        if self.ptm_file is not None:
            ptm = {"mode": "file", "file": self.ptm_file}
        else:
            ptm = {
                "mode": "inline",
                "states": self.ptm_states or [{"name": "WT", "modifications": []}],
            }
        return {
            "project": {"name": self.name, "outdir": str(self.outdir)},
            "input": {
                "protein": self.protein.to_dict(),
                "ptm": ptm,
                "mutations": {"mode": "inline", "states": self.mutation_states}
                if self.mutation_states
                else {"mode": "none"},
                "cleavage": self.cleavage,
            },
            "protocol": self.protocol,
            "analysis": self.analysis,
            "report": self.report,
            "sweep": self.sweep,
        }

    def write_yaml(self, path: str | Path) -> Path:
        """Write the concise user config YAML."""

        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(self.to_config(), sort_keys=False), encoding="utf-8")
        return output

    def compile(self, *, force: bool = True) -> LockedConfig:
        """Compile the recipe into project lock files."""

        normalized = normalize_config(self.to_config())
        resolved = resolve_presets(normalized)
        locked = compile_config(resolved)
        self.write_lock(locked, force=force)
        return locked

    def write_lock(self, locked: LockedConfig, *, force: bool = True) -> None:
        """Write a compiled lock produced from this experiment."""

        from idrptm.configuration import write_locked_config

        write_locked_config(locked, force=force)
