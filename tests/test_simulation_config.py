from __future__ import annotations

import pytest
from pydantic import ValidationError

from idrptm.schema import SimulationConfig


def test_save_every_steps_and_n_frames_derives_total_time() -> None:
    simulation = SimulationConfig(save_every_steps=7000, n_frames=1010)

    assert simulation.total_steps == 7_070_000
    assert simulation.frame_interval_ns == pytest.approx(0.07)
    assert simulation.total_time_ns == pytest.approx(70.7)


def test_total_time_and_frame_interval_derives_steps_and_frames() -> None:
    simulation = SimulationConfig(total_time_ns=10.0, frame_interval_ns=0.05)

    assert simulation.save_every_steps == 5000
    assert simulation.n_frames == 200
    assert simulation.total_steps == 1_000_000
    assert simulation.total_time_ns == pytest.approx(10.0)


def test_custom_dt_raises_validation_error() -> None:
    with pytest.raises(ValidationError, match="custom dt/integrator switching"):
        SimulationConfig(dt_ps=0.02)


def test_custom_integrator_raises_validation_error() -> None:
    with pytest.raises(ValidationError, match="custom dt/integrator switching"):
        SimulationConfig(integrator="custom")
