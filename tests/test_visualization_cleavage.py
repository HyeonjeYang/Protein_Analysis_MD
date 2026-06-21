from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from idrptm.visualization.cleavage import (
    event_schedule_figure,
    fragment_architecture_figure,
    plot_cut_number_trend,
    plot_fragment_length_distribution,
)


def _fragments() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "fragment_id": ["f1", "f2"],
            "original_start": [1, 5],
            "original_end": [4, 8],
        }
    )


def test_cleavage_architecture_and_lengths_plot() -> None:
    figs = [
        fragment_architecture_figure(_fragments()),
        plot_fragment_length_distribution(_fragments()),
    ]

    for fig in figs:
        assert fig.axes
        plt.close(fig)


def test_event_schedule_and_cut_number_trend_are_raw_points() -> None:
    events = pd.DataFrame(
        {"event_time_ns": [1.0, 2.5], "cut_number": [1, 2], "cut_after": [4, 8]}
    )
    fig = event_schedule_figure(events)
    assert list(fig.axes[0].lines[0].get_xdata()) == [1.0, 2.5]
    plt.close(fig)

    trend = pd.DataFrame({"cut_number": [0, 1, 2], "mean_rg": [1.0, 1.2, 1.4]})
    fig = plot_cut_number_trend(trend, "mean_rg", ylabel="Rg (nm)", title="Cut trend")
    assert fig.axes[0].collections
    plt.close(fig)
