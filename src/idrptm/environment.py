"""Environment metadata models and helpers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Environment:
    """Simulation environment settings."""

    pH: float = 7.4
    ionic_M: float = 0.15
    temperature_K: float = 298.0
    stickiness_scale: float = 1.0

    def metadata(self) -> dict[str, float | str]:
        """Return serializable metadata with interpretation hints."""

        return {
            "pH": self.pH,
            "ionic_M": self.ionic_M,
            "temperature_K": self.temperature_K,
            "stickiness_scale": self.stickiness_scale,
            "stickiness_note": (
                "non-default stickiness scaling is a model perturbation"
                if self.stickiness_scale != 1.0
                else "none"
            ),
        }
