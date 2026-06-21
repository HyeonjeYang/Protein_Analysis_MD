"""Residue-parameter writing for CALVADOS-compatible run directories."""

from __future__ import annotations

import csv
import importlib.util
import os
import sys
from dataclasses import dataclass
from pathlib import Path

RESIDUE_SOURCE_ENV_VARS = ("IDRPTM_CALVADOS_RESIDUES", "CALVADOS_RESIDUES_CSV")
REQUIRED_RESIDUE_COLUMNS = ("one", "three", "MW", "lambdas", "sigmas", "q", "bondlength")


@dataclass(frozen=True)
class ResidueParameterSet:
    """Reference to a CALVADOS-compatible residue parameter table."""

    name: str
    path: Path | None = None
    supports_ptms: tuple[str, ...] = ()


@dataclass(frozen=True)
class PhosphorylatedResidue:
    """Parameters needed to add or override one phosphorylated residue row."""

    ptm: str
    one: str
    three: str
    mw: float
    lambdas: float
    sigmas: float
    pka: float
    bondlength: float = 0.38

    def charge(self, ph: float) -> float:
        """Return the pH-dependent net charge used in CALVADOS pIDR examples."""

        return -1.0 - 1.0 / (1.0 + 10.0 ** (self.pka - ph))

    def row(self, ph: float, fieldnames: tuple[str, ...]) -> dict[str, str]:
        """Return a CSV row, preserving any extra source columns as blanks."""

        values = {fieldname: "" for fieldname in fieldnames}
        values.update(
            {
                "one": self.one,
                "three": self.three,
                "MW": _format_float(self.mw),
                "lambdas": _format_float(self.lambdas),
                "sigmas": _format_float(self.sigmas),
                "q": _format_float(self.charge(ph)),
                "bondlength": _format_float(self.bondlength),
            }
        )
        return values


@dataclass(frozen=True)
class ResidueWriteResult:
    """Result of writing a per-run residue parameter file."""

    source_file: Path
    output_file: Path
    ph: float
    charges: dict[str, dict[str, float | str]]

    def metadata(self) -> dict[str, object]:
        """Return JSON-serializable residue-parameter metadata."""

        return {
            "source_file": str(self.source_file),
            "output_file": str(self.output_file),
            "ph": self.ph,
            "ptm_charges": self.charges,
        }


PHOSPHORYLATED_RESIDUES = (
    PhosphorylatedResidue(
        ptm="pSer",
        one="B",
        three="SEP",
        mw=165.04,
        lambdas=0.0925,
        sigmas=0.601,
        pka=6.01,
    ),
    PhosphorylatedResidue(
        ptm="pThr",
        one="O",
        three="TPO",
        mw=179.07,
        lambdas=0.0013,
        sigmas=0.635,
        pka=6.30,
    ),
)


def default_parameter_sets() -> list[ResidueParameterSet]:
    """Return known parameter-set placeholders.

    Concrete files are supplied by the user or by an environment variable so
    upstream CALVADOS files are never mutated.
    """

    return [
        ResidueParameterSet(name="CALVADOS2", supports_ptms=()),
        ResidueParameterSet(name="pCALVADOS2", supports_ptms=("pSer", "pThr")),
    ]


def resolve_residue_source(path: str | Path | None = None) -> Path:
    """Resolve the base CALVADOS residue CSV from config or environment."""

    if path is not None:
        return _existing_file(Path(path), "configured residue parameter file")

    for variable in RESIDUE_SOURCE_ENV_VARS:
        value = os.environ.get(variable)
        if value:
            return _existing_file(Path(value), f"${variable}")

    installed_source = _installed_calvados_residue_csv()
    if installed_source is not None:
        return installed_source

    variables = ", ".join(RESIDUE_SOURCE_ENV_VARS)
    raise ValueError(
        "No base CALVADOS residue CSV was provided. Set calvados.residue_parameters "
        f"or one of: {variables}. If CALVADOS is installed, expected "
        "calvados/data/residues.csv to be importable."
    )


def phosphorylated_charge(ptm: str, ph: float) -> float:
    """Return the pH-dependent charge for pSer or pThr."""

    residue = _phosphorylated_by_ptm(ptm)
    return residue.charge(ph)


def write_residue_parameters(
    *,
    source_csv: str | Path | None,
    output_csv: str | Path,
    ph: float,
) -> ResidueWriteResult:
    """Write a per-run CALVADOS residues.csv with pSer/pThr rows updated."""

    source_file = resolve_residue_source(source_csv)
    output_file = Path(output_csv)
    rows, fieldnames = _read_residue_rows(source_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    rows_by_one = {row["one"]: index for index, row in enumerate(rows)}
    charges: dict[str, dict[str, float | str]] = {}
    for residue in PHOSPHORYLATED_RESIDUES:
        row = residue.row(ph, fieldnames)
        if residue.one in rows_by_one:
            rows[rows_by_one[residue.one]] = row
        else:
            rows_by_one[residue.one] = len(rows)
            rows.append(row)
        charges[residue.ptm] = {
            "one": residue.one,
            "three": residue.three,
            "pka": residue.pka,
            "charge": residue.charge(ph),
        }

    with output_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    return ResidueWriteResult(
        source_file=source_file,
        output_file=output_file,
        ph=ph,
        charges=charges,
    )


def _read_residue_rows(source_file: Path) -> tuple[list[dict[str, str]], tuple[str, ...]]:
    with source_file.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Residue CSV {source_file} has no header.")
        fieldnames = tuple(reader.fieldnames)
        missing = [column for column in REQUIRED_RESIDUE_COLUMNS if column not in fieldnames]
        if missing:
            missing_text = ", ".join(missing)
            raise ValueError(f"Residue CSV {source_file} is missing column(s): {missing_text}.")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Residue CSV {source_file} contains no residue rows.")
    return rows, fieldnames


def _phosphorylated_by_ptm(ptm: str) -> PhosphorylatedResidue:
    for residue in PHOSPHORYLATED_RESIDUES:
        if residue.ptm == ptm:
            return residue
    supported = ", ".join(residue.ptm for residue in PHOSPHORYLATED_RESIDUES)
    raise ValueError(f"Unsupported phosphorylated residue {ptm!r}; expected one of {supported}.")


def _existing_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"{label} does not exist or is not a file: {resolved}")
    return resolved


def _installed_calvados_residue_csv() -> Path | None:
    for entry in sys.path:
        root = Path(entry or ".").resolve() / "calvados"
        candidate = root / "data" / "residues.csv"
        if candidate.is_file():
            return candidate.resolve()

    spec = importlib.util.find_spec("calvados")
    if spec is None:
        return None
    locations = spec.submodule_search_locations
    if locations:
        roots = [Path(location) for location in locations]
    elif spec.origin:
        roots = [Path(spec.origin).parent]
    else:
        return None
    for root in roots:
        candidate = root / "data" / "residues.csv"
        if candidate.is_file():
            return candidate.resolve()
    return None


def _format_float(value: float) -> str:
    return f"{value:.15g}"
