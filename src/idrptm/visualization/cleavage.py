"""Cleavage/proteolysis visualization helpers."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from idrptm.plotting.plots import plot_cleavage_map, plot_event_schedule, plot_fragment_architecture


def plot_fragment_length_distribution(fragments: pd.DataFrame) -> plt.Figure:
    """Plot fragment length distribution."""

    lengths = _fragment_lengths(fragments)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(lengths, bins=max(1, min(20, len(lengths))))
    ax.set_xlabel("Fragment length (residues)")
    ax.set_ylabel("Fragment count")
    ax.set_title("Fragment length distribution")
    return fig


def plot_cut_site_event_scatter(events: pd.DataFrame) -> plt.Figure:
    """Plot cut-after position against raw event time."""

    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.scatter(events["event_time_ns"], events["cut_after"], s=28)
    ax.set_xlabel("Event time (ns)")
    ax.set_ylabel("Cut-after site (residue)")
    ax.set_title("Cleavage events (raw schedule)")
    return fig


def plot_cut_number_trend(
    table: pd.DataFrame,
    y: str,
    *,
    ylabel: str,
    title: str,
    show_visual_trend: bool = False,
) -> plt.Figure:
    """Plot raw values versus cut number; optional line is only a visual guide."""

    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.scatter(table["cut_number"], table[y], label="raw state values")
    if show_visual_trend:
        ordered = table.sort_values("cut_number")
        ax.plot(
            ordered["cut_number"],
            ordered[y],
            color="0.35",
            linewidth=1.5,
            label="visual guide",
        )
    ax.set_xlabel("Cut number (count)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} (quasi-dynamic cleavage states)")
    ax.legend(frameon=False)
    return fig


def plot_cleavage_architecture(
    sequence_length: int,
    cleavage_sites: pd.DataFrame,
    fragments: pd.DataFrame,
) -> plt.Figure:
    """Plot cleavage map and fragment architecture in one panel."""

    fig, axes = plt.subplots(2, 1, figsize=(8, 4.2), sharex=True)
    axes[0].hlines(0, 1, sequence_length, color="0.75", linewidth=6)
    if not cleavage_sites.empty:
        axes[0].vlines(cleavage_sites["cut_after"], -0.35, 0.35, color="tab:red", linewidth=2)
    for row_index, row in fragments.reset_index(drop=True).iterrows():
        axes[1].hlines(
            row_index,
            int(row["original_start"]),
            int(row["original_end"]),
            linewidth=5,
            color="tab:blue",
        )
    axes[1].set_yticks(range(len(fragments)), fragments["fragment_id"].astype(str).tolist())
    axes[0].set_yticks([])
    axes[0].set_ylim(-0.8, 0.8)
    axes[0].set_title("Cleavage map")
    axes[1].set_title("Fragment architecture")
    axes[1].set_xlabel("Original residue position (residue)")
    axes[1].set_xlim(1, max(sequence_length, 1))
    return fig


def cleavage_map_figure(sequence_length: int, cleavage_sites: pd.DataFrame) -> plt.Figure:
    """Wrapper for public visualization namespace."""

    return plot_cleavage_map(sequence_length, cleavage_sites)


def fragment_architecture_figure(fragments: pd.DataFrame) -> plt.Figure:
    """Wrapper for public visualization namespace."""

    return plot_fragment_architecture(fragments)


def event_schedule_figure(events: pd.DataFrame) -> plt.Figure:
    """Return raw cleavage event schedule plot."""

    return plot_event_schedule(events)


def _fragment_lengths(fragments: pd.DataFrame) -> list[int]:
    if {"original_start", "original_end"}.issubset(fragments.columns):
        return (fragments["original_end"] - fragments["original_start"] + 1).astype(int).tolist()
    if "sequence" in fragments:
        return fragments["sequence"].astype(str).str.len().astype(int).tolist()
    return []
