"""Strict visualization smoothing policy."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SmoothingRule:
    """Policy for one data product."""

    target: str
    default_enabled: bool
    allowed: bool
    warning_required: bool
    visualization_only: bool
    raw_is_quantitative_default: bool
    note: str


_RULES: dict[str, SmoothingRule] = {
    "ps": SmoothingRule("ps", True, True, False, True, True, "Log-space trend display is allowed."),
    "rs": SmoothingRule("rs", True, True, False, True, True, "Log-space trend display is allowed."),
    "energy": SmoothingRule(
        "energy",
        True,
        True,
        False,
        True,
        True,
        "Rolling visual drift checks are allowed.",
    ),
    "rg": SmoothingRule("rg", False, True, False, True, True, "Rolling visual trend only."),
    "ree": SmoothingRule("ree", False, True, False, True, True, "Rolling visual trend only."),
    "density_profile": SmoothingRule(
        "density_profile",
        True,
        True,
        False,
        True,
        True,
        "Display smoothing is allowed; concentration tables remain raw/default-profile based.",
    ),
    "contact_map": SmoothingRule(
        "contact_map",
        False,
        True,
        True,
        True,
        True,
        "Matrix smoothing is visualization-only and disabled by default.",
    ),
    "delta_contact_map": SmoothingRule(
        "delta_contact_map",
        False,
        True,
        True,
        True,
        True,
        "Smooth only after raw delta calculation, for display only.",
    ),
    "contact_lifetime": SmoothingRule(
        "contact_lifetime",
        False,
        False,
        True,
        False,
        True,
        "Do not smooth lifetime calculations.",
    ),
    "binary_contact_time_series": SmoothingRule(
        "binary_contact_time_series",
        False,
        False,
        True,
        False,
        True,
        "Do not smooth binary state series.",
    ),
    "event_schedule": SmoothingRule(
        "event_schedule",
        False,
        False,
        True,
        False,
        True,
        "Do not smooth cleavage/Poisson event schedules.",
    ),
    "cut_number": SmoothingRule(
        "cut_number",
        False,
        False,
        True,
        False,
        True,
        "Do not smooth discrete cut numbers.",
    ),
    "sequence_features": SmoothingRule(
        "sequence_features",
        False,
        False,
        True,
        False,
        True,
        "Do not smooth sequence-derived features.",
    ),
}


def smoothing_rule(target: str) -> SmoothingRule:
    """Return the policy rule for a target, defaulting to conservative off."""

    key = target.lower()
    return _RULES.get(
        key,
        SmoothingRule(
            key,
            False,
            True,
            True,
            True,
            True,
            (
                "Unknown target: smoothing is optional, visualization-only, "
                "and raw metrics remain default."
            ),
        ),
    )


def validate_smoothing_request(target: str, config: dict[str, object] | None) -> dict[str, object]:
    """Validate a smoothing block and return metadata suitable for summaries."""

    rule = smoothing_rule(target)
    payload = dict(config or {})
    enabled = bool(payload.get("enabled", rule.default_enabled))
    if enabled and not rule.allowed:
        raise ValueError(f"Smoothing is not allowed for {target}: {rule.note}")
    if enabled and rule.visualization_only:
        payload["visualization_only"] = True
    payload.setdefault("enabled", enabled)
    payload["raw_is_quantitative_default"] = rule.raw_is_quantitative_default
    payload["warning_required"] = bool(enabled and rule.warning_required)
    payload["policy_note"] = rule.note
    return payload


def conservative_smoothing_defaults() -> dict[str, dict[str, object]]:
    """Return the conservative default smoothing policy block."""

    return {
        "policy": {"name": "conservative", "raw_is_source_of_truth": True},
        "ps": {
            "enabled": True,
            "method": "logspace",
            "window_log10": 0.2,
            "min_points": 5,
            "robust": True,
            "use_for_fit": False,
        },
        "rs": {
            "enabled": True,
            "method": "logspace",
            "window_log10": 0.2,
            "min_points": 5,
            "robust": True,
            "use_for_fit": False,
        },
        "energy": {"enabled": True, "method": "rolling", "window": 25, "visualization_only": True},
        "rg": {"enabled": False, "method": "rolling", "window": 25, "visualization_only": True},
        "ree": {"enabled": False, "method": "rolling", "window": 25, "visualization_only": True},
        "density_profile": {
            "enabled": True,
            "method": "gaussian",
            "sigma_bins": 1.0,
            "visualization_only": True,
        },
        "contact_map": {
            "enabled": False,
            "method": "gaussian",
            "sigma": 1.0,
            "visualization_only": True,
        },
        "delta_contact_map": {
            "enabled": False,
            "method": "gaussian",
            "sigma": 1.0,
            "visualization_only": True,
        },
    }


def event_schedules_are_never_smoothed() -> bool:
    """Tiny explicit predicate used by tests and report code."""

    return not smoothing_rule("event_schedule").allowed
