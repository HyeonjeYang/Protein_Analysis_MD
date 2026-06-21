"""Pre-simulation sequence cleavage and proteolysis perturbations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from idrptm.ptm import AppliedPTM
from idrptm.schema import CleavageProduct, CleavageSet, ProteaseRule
from idrptm.sequence import FASTA_WRAP, normalize_raw_sequence, sanitize_identifier


@dataclass(frozen=True)
class CleavageSite:
    """One selected cleavage site using one-based residue numbering."""

    cut_after: int
    n_terminal_residue: str
    c_terminal_residue: str
    rule: str


@dataclass(frozen=True)
class CleavageState:
    """One pre-cleaved sequence state."""

    name: str
    cut_number: int
    cuts: tuple[int, ...]
    sites: tuple[CleavageSite, ...]
    products: tuple[CleavageProduct, ...]
    event_time_ns: float | None = None


def protease_candidate_cuts(sequence: str, rule: str | ProteaseRule) -> tuple[int, ...]:
    """Return one-based cut positions for a built-in protease rule."""

    normalized = normalize_raw_sequence(sequence)
    rule_name = _rule_name(rule)
    if rule_name == "trypsin_simple":
        return _cuts_after_residues(normalized, {"K", "R"}, block_before={"P"})
    if rule_name == "lysc_simple":
        return _cuts_after_residues(normalized, {"K"})
    if rule_name == "argc_simple":
        return _cuts_after_residues(normalized, {"R"})
    if rule_name == "chymotrypsin_high_simple":
        return _cuts_after_residues(normalized, {"F", "Y", "W"}, block_before={"P"})
    if rule_name == "chymotrypsin_low_simple":
        return _cuts_after_residues(normalized, {"F", "Y", "W", "L", "M", "H"}, block_before={"P"})
    if rule_name == "cnbr_simple":
        return _cuts_after_residues(normalized, {"M"})
    if rule_name == "glu_c_simple":
        return _cuts_after_residues(normalized, {"E"})
    if rule_name == "asp_n_simple":
        return _cuts_before_residues(normalized, {"D"})
    if rule_name == "tev_simple":
        return _tev_cuts(normalized)
    raise ValueError(f"Unsupported protease rule: {rule_name}")


def generate_cleavage_states(
    sequence: str,
    cleavage_set: CleavageSet,
    *,
    ptm_sites: tuple[AppliedPTM, ...] = (),
) -> tuple[CleavageState, ...]:
    """Generate pre-cleaved sequence states from a cleavage set."""

    normalized = normalize_raw_sequence(sequence)
    rng = np.random.default_rng(cleavage_set.seed)
    if cleavage_set.mode == "none":
        return (
            _state_from_cuts(
                normalized,
                cleavage_set,
                cuts=(),
                cut_number=0,
                ptm_sites=ptm_sites,
            ),
        )
    if cleavage_set.mode == "poisson":
        events = poisson_cleavage_events(normalized, cleavage_set, rng=rng)
        return tuple(
            _state_from_cuts(
                normalized,
                cleavage_set,
                cuts=tuple(sorted(event.cut_after for event in events[:index])),
                cut_number=index,
                ptm_sites=ptm_sites,
                event_time_ns=events[index - 1].event_time_ns,
            )
            for index in range(1, len(events) + 1)
        )
    cuts = _select_cuts(normalized, cleavage_set, rng)
    if cleavage_set.mode == "sequential" or cleavage_set.sequential_series:
        ordered = _ordered_cuts(cuts, cleavage_set.order, rng)
        if cleavage_set.n_cuts is not None:
            ordered = ordered[: cleavage_set.n_cuts]
        return tuple(
            _state_from_cuts(
                normalized,
                cleavage_set,
                cuts=tuple(sorted(ordered[:index])),
                cut_number=index,
                ptm_sites=ptm_sites,
            )
            for index in range(1, len(ordered) + 1)
        )
    return (
        _state_from_cuts(
            normalized,
            cleavage_set,
            cuts=cuts,
            cut_number=len(cuts),
            ptm_sites=ptm_sites,
        ),
    )


def fragments_from_cuts(
    sequence: str,
    cuts: tuple[int, ...] | list[int],
    *,
    missed_cleavages: int = 0,
    min_fragment_length: int = 1,
    preserve_ptm_mapping: bool = True,
    ptm_sites: tuple[AppliedPTM, ...] = (),
    prefix: str = "fragment",
) -> tuple[CleavageProduct, ...]:
    """Return fragment products for selected cuts."""

    normalized = normalize_raw_sequence(sequence)
    boundaries = [0, *sorted(_valid_cut_set(cuts, len(normalized))), len(normalized)]
    base_ranges = list(zip(boundaries[:-1], boundaries[1:], strict=True))
    products: list[CleavageProduct] = []
    fragment_index = 1
    for start_index in range(len(base_ranges)):
        for end_index in range(
            start_index,
            min(start_index + missed_cleavages + 1, len(base_ranges)),
        ):
            start = base_ranges[start_index][0]
            end = base_ranges[end_index][1]
            if end - start < min_fragment_length:
                continue
            original_indices = list(range(start + 1, end + 1))
            products.append(
                CleavageProduct(
                    fragment_id=f"{prefix}_{fragment_index}",
                    sequence=normalized[start:end],
                    original_start=start + 1,
                    original_end=end,
                    original_indices=original_indices,
                    ptm_sites_1based=_fragment_ptm_sites(
                        original_indices,
                        ptm_sites,
                        preserve=preserve_ptm_mapping,
                    ),
                )
            )
            fragment_index += 1
    return tuple(products)


def write_fragments_fasta(
    products: tuple[CleavageProduct, ...] | list[CleavageProduct],
    path: str | Path,
) -> Path:
    """Write cleavage products to a FASTA file."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for product in products:
        records.append(_format_fragment_fasta(product).rstrip())
    output_path.write_text("\n".join(records) + "\n", encoding="utf-8")
    return output_path


@dataclass(frozen=True)
class CleavageEvent:
    """One scheduled quasi-dynamic cleavage event."""

    event_index: int
    event_time_ns: float
    cut_after: int


def poisson_cleavage_events(
    sequence: str,
    cleavage_set: CleavageSet,
    *,
    rng: np.random.Generator | None = None,
) -> tuple[CleavageEvent, ...]:
    """Generate a reproducible Poisson cleavage event schedule."""

    generator = rng or np.random.default_rng(cleavage_set.seed)
    normalized = normalize_raw_sequence(sequence)
    candidates = _poisson_candidates(normalized, cleavage_set)
    max_time = cleavage_set.max_time_ns or 0.0
    max_cuts = cleavage_set.max_cuts if cleavage_set.max_cuts is not None else cleavage_set.n_cuts
    max_cuts = max_cuts if max_cuts is not None else len(candidates)
    if max_time <= 0 or max_cuts <= 0:
        return ()
    site_rate = _site_rate(cleavage_set, len(candidates))
    current_time = 0.0
    selected: list[int] = []
    events: list[CleavageEvent] = []
    remaining = list(candidates)
    while remaining and len(events) < max_cuts:
        total_rate = site_rate * len(remaining)
        current_time += float(generator.exponential(1.0 / total_rate))
        if current_time > max_time:
            break
        cut = int(generator.choice(remaining))
        trial = sorted([*selected, cut])
        if _base_fragments_satisfy_min_length(
            trial,
            sequence_length=len(normalized),
            min_fragment_length=cleavage_set.min_fragment_length,
        ):
            selected.append(cut)
            events.append(
                CleavageEvent(
                    event_index=len(events) + 1,
                    event_time_ns=current_time,
                    cut_after=cut,
                )
            )
        remaining = [candidate for candidate in remaining if candidate not in selected]
    return tuple(events)


def _format_fragment_fasta(product: CleavageProduct) -> str:
    header = (
        f"{sanitize_identifier(product.fragment_id)} "
        f"original_range={product.original_start}-{product.original_end}"
    )
    sequence = normalize_raw_sequence(product.sequence)
    lines = [f">{header}"]
    lines.extend(
        sequence[index : index + FASTA_WRAP]
        for index in range(0, len(sequence), FASTA_WRAP)
    )
    return "\n".join(lines) + "\n"


def _select_cuts(
    sequence: str,
    cleavage_set: CleavageSet,
    rng: np.random.Generator,
) -> tuple[int, ...]:
    if cleavage_set.mode in {"protease", "enzyme"}:
        assert cleavage_set.protease is not None
        candidates = protease_candidate_cuts(sequence, cleavage_set.protease)
    elif cleavage_set.mode == "manual":
        candidates = tuple(cleavage_set.manual_cuts)
    elif cleavage_set.mode in {"random", "random_ncuts"}:
        candidates = tuple(range(1, len(sequence)))
    elif cleavage_set.mode == "sequential":
        candidates = _sequential_candidates(sequence, cleavage_set)
    elif cleavage_set.mode == "end_trimming":
        candidates = _end_trimming_cuts(sequence, cleavage_set)
    else:
        raise ValueError(f"Unsupported cleavage mode: {cleavage_set.mode}")

    candidates = _valid_cut_set(candidates, len(sequence))
    if cleavage_set.cleavage_probability < 1.0:
        candidates = tuple(
            cut for cut in candidates if rng.random() <= cleavage_set.cleavage_probability
        )
    if cleavage_set.mode in {"random", "random_ncuts"}:
        n_cuts = cleavage_set.n_cuts if cleavage_set.n_cuts is not None else 1
        return _random_cuts(
            candidates,
            n_cuts=n_cuts,
            sequence_length=len(sequence),
            min_fragment_length=cleavage_set.min_fragment_length,
            rng=rng,
        )
    if cleavage_set.n_cuts is not None:
        candidates = _ordered_cuts(candidates, cleavage_set.order, rng)[: cleavage_set.n_cuts]
    return tuple(sorted(candidates))


def _sequential_candidates(sequence: str, cleavage_set: CleavageSet) -> tuple[int, ...]:
    if cleavage_set.protease is not None:
        return protease_candidate_cuts(sequence, cleavage_set.protease)
    if cleavage_set.manual_cuts:
        return tuple(cleavage_set.manual_cuts)
    return tuple(range(1, len(sequence)))


def _state_from_cuts(
    sequence: str,
    cleavage_set: CleavageSet,
    *,
    cuts: tuple[int, ...],
    cut_number: int,
    ptm_sites: tuple[AppliedPTM, ...],
    event_time_ns: float | None = None,
) -> CleavageState:
    products = fragments_from_cuts(
        sequence,
        cuts,
        missed_cleavages=cleavage_set.missed_cleavages,
        min_fragment_length=cleavage_set.min_fragment_length,
        preserve_ptm_mapping=cleavage_set.preserve_ptm_mapping,
        ptm_sites=ptm_sites,
        prefix=f"{cleavage_set.name}_cut{cut_number}_fragment",
    )
    return CleavageState(
        name=f"{cleavage_set.name}_cut{cut_number}",
        cut_number=cut_number,
        cuts=tuple(sorted(cuts)),
        sites=_sites_from_cuts(sequence, cuts, _rule_name(cleavage_set.protease)),
        products=products,
        event_time_ns=event_time_ns,
    )


def _cuts_after_residues(
    sequence: str,
    residues: set[str],
    *,
    block_before: set[str] | None = None,
) -> tuple[int, ...]:
    block = block_before or set()
    cuts = []
    for index, residue in enumerate(sequence[:-1]):
        if residue in residues and sequence[index + 1] not in block:
            cuts.append(index + 1)
    return tuple(cuts)


def _cuts_before_residues(sequence: str, residues: set[str]) -> tuple[int, ...]:
    return tuple(
        index
        for index, residue in enumerate(sequence[1:], start=1)
        if residue in residues
    )


def _tev_cuts(sequence: str) -> tuple[int, ...]:
    cuts = []
    for index in range(0, len(sequence) - 6):
        if (
            sequence[index] == "E"
            and sequence[index + 3] == "Y"
            and sequence[index + 5] == "Q"
            and sequence[index + 6] in {"G", "S"}
        ):
            cuts.append(index + 6)
    return tuple(cuts)


def _end_trimming_cuts(sequence: str, cleavage_set: CleavageSet) -> tuple[int, ...]:
    assert cleavage_set.terminus is not None
    assert cleavage_set.step_size is not None
    max_removed = (
        cleavage_set.max_removed
        if cleavage_set.max_removed is not None
        else len(sequence) - 1
    )
    values = range(
        cleavage_set.step_size,
        min(max_removed, len(sequence) - 1) + 1,
        cleavage_set.step_size,
    )
    if cleavage_set.terminus == "N":
        return tuple(values)
    return tuple(len(sequence) - value for value in values)


def _poisson_candidates(sequence: str, cleavage_set: CleavageSet) -> tuple[int, ...]:
    if cleavage_set.candidate_sites == "enzyme" and cleavage_set.protease is not None:
        return protease_candidate_cuts(sequence, cleavage_set.protease)
    if isinstance(cleavage_set.candidate_sites, list):
        return _valid_cut_set(cleavage_set.candidate_sites, len(sequence))
    if cleavage_set.manual_cuts:
        return _valid_cut_set(cleavage_set.manual_cuts, len(sequence))
    return tuple(range(1, len(sequence)))


def _site_rate(cleavage_set: CleavageSet, n_candidates: int) -> float:
    if cleavage_set.site_rate_per_ns is not None:
        return cleavage_set.site_rate_per_ns
    if cleavage_set.global_rate_per_ns is not None and n_candidates > 0:
        return cleavage_set.global_rate_per_ns / n_candidates
    return 0.01


def _valid_cut_set(cuts: tuple[int, ...] | list[int], sequence_length: int) -> tuple[int, ...]:
    valid = sorted({int(cut) for cut in cuts if 0 < int(cut) < sequence_length})
    return tuple(valid)


def _ordered_cuts(
    cuts: tuple[int, ...],
    order: str,
    rng: np.random.Generator,
) -> tuple[int, ...]:
    if order == "n_to_c":
        return tuple(sorted(cuts))
    if order == "c_to_n":
        return tuple(sorted(cuts, reverse=True))
    if order == "random":
        shuffled = list(cuts)
        rng.shuffle(shuffled)
        return tuple(shuffled)
    raise ValueError(f"Unsupported cleavage order: {order}")


def _random_cuts(
    candidates: tuple[int, ...],
    *,
    n_cuts: int,
    sequence_length: int,
    min_fragment_length: int,
    rng: np.random.Generator,
) -> tuple[int, ...]:
    shuffled = list(candidates)
    rng.shuffle(shuffled)
    selected: list[int] = []
    for cut in shuffled:
        trial = sorted([*selected, cut])
        if _base_fragments_satisfy_min_length(
            trial,
            sequence_length=sequence_length,
            min_fragment_length=min_fragment_length,
        ):
            selected.append(cut)
        if len(selected) >= n_cuts:
            break
    return tuple(sorted(selected))


def _base_fragments_satisfy_min_length(
    cuts: list[int],
    *,
    sequence_length: int,
    min_fragment_length: int,
) -> bool:
    boundaries = [0, *sorted(cuts), sequence_length]
    return all(
        (end - start) >= min_fragment_length
        for start, end in zip(boundaries[:-1], boundaries[1:], strict=True)
    )


def _fragment_ptm_sites(
    original_indices: list[int],
    ptm_sites: tuple[AppliedPTM, ...],
    *,
    preserve: bool,
) -> list[int]:
    if not preserve:
        return []
    index_set = set(original_indices)
    return [
        site.biological_position
        for site in ptm_sites
        if site.biological_position in index_set
    ]


def _sites_from_cuts(
    sequence: str,
    cuts: tuple[int, ...] | list[int],
    rule: str,
) -> tuple[CleavageSite, ...]:
    return tuple(
        CleavageSite(
            cut_after=cut,
            n_terminal_residue=sequence[cut - 1],
            c_terminal_residue=sequence[cut],
            rule=rule,
        )
        for cut in sorted(cuts)
    )


def _rule_name(rule: str | ProteaseRule | None) -> str:
    if rule is None:
        return "manual"
    if isinstance(rule, ProteaseRule):
        return rule.name
    return rule
