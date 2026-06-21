"""UniProt/Swiss-Prot sequence search, fetch, and cache helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"


class UniProtError(RuntimeError):
    """Base error for UniProt lookup failures."""


class AmbiguousUniProtQuery(UniProtError):
    """Raised when a query has multiple plausible candidates."""


class UniProtNoMatch(UniProtError):
    """Raised when no candidate was found."""


@dataclass(frozen=True)
class UniProtCandidate:
    """Small ranked UniProt candidate record."""

    accession: str
    entry_name: str
    reviewed: bool
    protein_name: str
    gene_names: tuple[str, ...]
    organism: str
    length: int
    sequence: str | None = None
    score: float = 0.0

    def to_dict(self) -> dict[str, object]:
        return {
            "accession": self.accession,
            "entry_name": self.entry_name,
            "reviewed": self.reviewed,
            "protein_name": self.protein_name,
            "gene_names": list(self.gene_names),
            "organism": self.organism,
            "length": self.length,
            "sequence": self.sequence,
            "score": self.score,
        }


@dataclass(frozen=True)
class UniProtFetchResult:
    """Fetched sequence and metadata paths."""

    candidate: UniProtCandidate
    fasta_path: Path
    metadata_path: Path
    metadata: dict[str, object]


def search_uniprot(
    query: str,
    *,
    reviewed: bool = True,
    organism: str | None = None,
    size: int = 10,
    session: Any | None = None,
) -> list[UniProtCandidate]:
    """Search UniProtKB and return ranked candidates."""

    client = session or requests
    params = {
        "query": _build_query(query, reviewed=reviewed, organism=organism),
        "format": "json",
        "size": str(size),
        "fields": "accession,id,reviewed,protein_name,gene_names,organism_name,length,sequence",
    }
    response = client.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
    response.raise_for_status()
    payload = response.json()
    candidates = [
        _candidate_from_entry(entry, query, organism)
        for entry in payload.get("results", [])
    ]
    return sorted(candidates, key=lambda candidate: candidate.score, reverse=True)


def fetch_sequence(
    query: str,
    *,
    accession: str | None = None,
    reviewed: bool = True,
    organism: str | None = None,
    interactive: bool = False,
    yes: bool = False,
    refresh: bool = False,
    cache_dir: str | Path = "data/sequences",
    session: Any | None = None,
) -> UniProtFetchResult:
    """Fetch one UniProt sequence and cache FASTA plus metadata."""

    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    client = session or requests

    if accession:
        candidate = _fetch_accession(accession, query=query, organism=organism, session=client)
    else:
        candidates = search_uniprot(query, reviewed=reviewed, organism=organism, session=client)
        if not candidates:
            raise UniProtNoMatch(f"No UniProt candidates found for query {query!r}.")
        candidate = _select_candidate(candidates, query=query, interactive=interactive, yes=yes)

    fasta_path, metadata_path = _cache_paths(cache, candidate)
    if fasta_path.exists() and metadata_path.exists() and not refresh:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        return UniProtFetchResult(
            candidate=candidate,
            fasta_path=fasta_path,
            metadata_path=metadata_path,
            metadata=metadata,
        )

    metadata = _metadata(candidate, query=query)
    fasta_path.write_text(_format_fasta(candidate), encoding="utf-8")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return UniProtFetchResult(
        candidate=candidate,
        fasta_path=fasta_path,
        metadata_path=metadata_path,
        metadata=metadata,
    )


def format_candidates(candidates: list[UniProtCandidate]) -> str:
    """Return a concise candidate table."""

    lines = ["rank accession entry_name reviewed genes organism length protein_name"]
    for index, candidate in enumerate(candidates, start=1):
        genes = ",".join(candidate.gene_names) or "-"
        lines.append(
            f"{index} {candidate.accession} {candidate.entry_name} "
            f"{candidate.reviewed} {genes} {candidate.organism} "
            f"{candidate.length} {candidate.protein_name}"
        )
    return "\n".join(lines)


def _build_query(query: str, *, reviewed: bool, organism: str | None) -> str:
    parts = [f"({query})"]
    if reviewed:
        parts.append("(reviewed:true)")
    if organism:
        parts.append(f'(organism_name:"{organism}")')
    return " AND ".join(parts)


def _candidate_from_entry(
    entry: dict[str, Any],
    query: str,
    organism: str | None,
) -> UniProtCandidate:
    accession = entry.get("primaryAccession", "")
    entry_name = entry.get("uniProtkbId", "")
    protein_name = _recommended_name(entry)
    genes = tuple(gene.get("geneName", {}).get("value", "") for gene in entry.get("genes", []))
    genes = tuple(gene for gene in genes if gene)
    entry_organism = entry.get("organism", {}).get("scientificName", "")
    sequence_data = entry.get("sequence", {}) or {}
    reviewed = entry.get("entryType", "").lower().startswith("uniprotkb reviewed")
    candidate = UniProtCandidate(
        accession=accession,
        entry_name=entry_name,
        reviewed=reviewed,
        protein_name=protein_name,
        gene_names=genes,
        organism=entry_organism,
        length=int(sequence_data.get("length") or 0),
        sequence=sequence_data.get("value"),
    )
    return candidate.__class__(
        **{**candidate.to_dict(), "score": _score_candidate(candidate, query, organism)}
    )


def _fetch_accession(
    accession: str,
    *,
    query: str,
    organism: str | None,
    session: Any,
) -> UniProtCandidate:
    response = session.get(UNIPROT_ENTRY_URL.format(accession=accession), timeout=30)
    response.raise_for_status()
    return _candidate_from_entry(response.json(), query, organism)


def _recommended_name(entry: dict[str, Any]) -> str:
    description = entry.get("proteinDescription", {})
    recommended = description.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    if full_name.get("value"):
        return str(full_name["value"])
    submission = description.get("submissionNames", [])
    if submission:
        return str(submission[0].get("fullName", {}).get("value", ""))
    return ""


def _score_candidate(candidate: UniProtCandidate, query: str, organism: str | None) -> float:
    query_upper = query.upper()
    score = 0.0
    if candidate.accession.upper() == query_upper:
        score += 100
    if query_upper in {gene.upper() for gene in candidate.gene_names}:
        score += 70
    if query_upper in candidate.entry_name.upper():
        score += 30
    if query_upper in candidate.protein_name.upper():
        score += 25
    if candidate.reviewed:
        score += 20
    if organism and organism.upper() == candidate.organism.upper():
        score += 15
    if candidate.sequence:
        score += 5
    return score


def _select_candidate(
    candidates: list[UniProtCandidate],
    *,
    query: str,
    interactive: bool,
    yes: bool,
) -> UniProtCandidate:
    if len(candidates) == 1:
        return candidates[0]
    if yes and candidates[0].score > candidates[1].score + 30:
        return candidates[0]
    if interactive:
        print(format_candidates(candidates))
        choice = input("Select UniProt candidate rank: ").strip()
        try:
            index = int(choice) - 1
        except ValueError as exc:
            raise AmbiguousUniProtQuery("Selection must be a candidate rank.") from exc
        if index < 0 or index >= len(candidates):
            raise AmbiguousUniProtQuery("Selected rank is out of range.")
        return candidates[index]
    raise AmbiguousUniProtQuery(
        f"Query {query!r} is ambiguous. Use --interactive, --accession, or --yes "
        "when there is a strong unique candidate."
    )


def _cache_paths(cache: Path, candidate: UniProtCandidate) -> tuple[Path, Path]:
    stem = f"{candidate.accession}_{candidate.entry_name}".replace("/", "_")
    return cache / f"{stem}.fasta", cache / f"{stem}.metadata.json"


def _metadata(candidate: UniProtCandidate, *, query: str) -> dict[str, object]:
    return {
        **candidate.to_dict(),
        "retrieval_timestamp": datetime.now(timezone.utc).isoformat(),
        "query_used": query,
        "source": "UniProtKB/Swiss-Prot" if candidate.reviewed else "UniProtKB",
    }


def _format_fasta(candidate: UniProtCandidate) -> str:
    sequence = candidate.sequence or ""
    header = f">{candidate.accession}|{candidate.entry_name} {candidate.protein_name}"
    chunks = [sequence[index : index + 80] for index in range(0, len(sequence), 80)]
    return "\n".join([header, *chunks]) + "\n"
