"""Simple enzyme cleavage-rule registry for MVP proteolysis designs."""

from __future__ import annotations

from dataclasses import dataclass

from idrptm.cleavage import protease_candidate_cuts


@dataclass(frozen=True)
class EnzymeRule:
    """Rule-based cleavage enzyme description."""

    name: str
    description: str


BUILTIN_ENZYMES: dict[str, EnzymeRule] = {
    "trypsin_simple": EnzymeRule("trypsin_simple", "cleave after K/R unless next residue is P"),
    "lysc_simple": EnzymeRule("lysc_simple", "cleave after K"),
    "argc_simple": EnzymeRule("argc_simple", "cleave after R"),
    "chymotrypsin_high_simple": EnzymeRule(
        "chymotrypsin_high_simple",
        "cleave after F/Y/W unless next residue is P",
    ),
    "chymotrypsin_low_simple": EnzymeRule(
        "chymotrypsin_low_simple",
        "planned low-specificity chymotrypsin rule",
    ),
    "cnbr_simple": EnzymeRule("cnbr_simple", "cleave after M"),
    "glu_c_simple": EnzymeRule("glu_c_simple", "planned Glu-C simple rule"),
    "asp_n_simple": EnzymeRule("asp_n_simple", "planned Asp-N simple rule"),
    "tev_simple": EnzymeRule("tev_simple", "E-X-X-Y-X-Q-[G/S], cleave after Q"),
}


def cleavage_sites(sequence: str, enzyme: str) -> tuple[int, ...]:
    """Return cleavage sites for a built-in enzyme."""

    return protease_candidate_cuts(sequence, enzyme)
