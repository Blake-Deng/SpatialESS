#!/usr/bin/env python3
"""Generate the two publication figures from bundled ovarian source tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmark" / "ovarian_20260904" / "source_data"
OUT = ROOT / "figures"

OFFICIAL = "#D55E00"
ESS = "#0072B2"
SPARKLE = "#009E73"
RAW = "#4D4D4D"
LOST = "#CC6677"
GAINED = "#4477AA"
SHARED = "#228833"
LIGHT = "#E8E8E8"
TEXT = "#202124"


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.edgecolor": "#555555",
            "axes.linewidth": 0.7,
            "text.color": TEXT,
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.07,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def clean_axes(ax, grid: bool = True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.55, alpha=0.75)
        ax.set_axisbelow(True)


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=450, bbox_inches="tight")
    plt.close(fig)


def figure_official_vs_ess() -> None:
    resources = pd.read_csv(DATA / "resource_summary.tsv", sep="\t")
    records = pd.read_csv(
        DATA / "official_vs_spatialess_records.tsv.gz", sep="\t"
    )
    exact = pd.read_csv(DATA / "official_vs_spatialess.tsv", sep="\t").iloc[0]

    fig = plt.figure(figsize=(7.2, 5.8), facecolor="white")
    grid = fig.add_gridspec(
        2, 2, height_ratios=[1.02, 1], width_ratios=[1.03, 1],
        hspace=0.38, wspace=0.33
    )

    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1)
    ax.axis("off")
    panel_label(ax, "A")
    ax.set_title("Matched benchmark scale", pad=8)
    summary = [
        ("16,247", "cells"),
        ("967", "signaling genes"),
        ("1,828", "ligand-receptor pairs"),
        ("100", "permutations"),
    ]
    y_positions = [0.82, 0.62, 0.42, 0.22]
    for (value, label), ypos in zip(summary, y_positions):
        ax.text(
            0.07, ypos, value, color=ESS, fontsize=13,
            fontweight="bold", va="center", ha="left"
        )
        ax.text(0.43, ypos, label, fontsize=8.2, va="center", ha="left")
        ax.plot([0.07, 0.93], [ypos - 0.09, ypos - 0.09], color=LIGHT, lw=0.8)
    ax.text(
        0.07, 0.02, "Same cells, labels, coordinates, LR set, seed and parameters",
        fontsize=6.9, color="#555555", va="bottom"
    )

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    ax.set_title("End-to-end core runtime", pad=8)
    order = ["Official SpatialCellChat V3", "SpatialESS"]
    values = resources.set_index("method").loc[order, "core_seconds"]
    y = np.arange(2)
    ax.hlines(y, 10, values, color=[OFFICIAL, ESS], linewidth=4.5, alpha=0.7)
    ax.scatter(values, y, s=48, color=[OFFICIAL, ESS], zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(10, 12000)
    ax.set_yticks(y, ["Official V3", "SpatialESS"])
    ax.invert_yaxis()
    ax.set_xlabel("Runtime (seconds, log scale)")
    for yi, value in zip(y, values):
        label = "1 h 49 min" if value > 1000 else f"{value:.1f} s"
        ax.text(value * 1.12, yi, label, va="center", fontsize=7.5)
    ax.text(
        0.98, 0.08, f"{values.iloc[0] / values.iloc[1]:.1f}x faster",
        transform=ax.transAxes, ha="right", color=ESS, fontweight="bold"
    )
    clean_axes(ax, grid=False)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.55, alpha=0.75)

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    ax.set_title("Peak resident memory", pad=8)
    memory = resources.set_index("method").loc[order, "peak_rss_gib"]
    bars = ax.bar(
        ["Official V3", "SpatialESS"], memory,
        color=[OFFICIAL, ESS], width=0.58, edgecolor="white", linewidth=0.8
    )
    ax.set_ylabel("Peak RSS (GiB)")
    ax.set_ylim(0, max(memory) * 1.28)
    for bar, value in zip(bars, memory):
        ax.text(
            bar.get_x() + bar.get_width() / 2, value + 0.10,
            f"{value:.2f}", ha="center", va="bottom", fontweight="bold"
        )
    reduction = 100 * (1 - memory.iloc[1] / memory.iloc[0])
    ax.text(
        0.98, 0.91, f"{reduction:.1f}% lower",
        transform=ax.transAxes, ha="right", color=ESS, fontweight="bold"
    )
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    ax.set_title("Communication probability concordance", pad=8)
    x = records["official_probability"].to_numpy()
    yv = records["spatialess_probability"].to_numpy()
    positive = np.concatenate([x[x > 0], yv[yv > 0]])
    floor = positive.min() * 0.7
    xplot = np.maximum(x, floor)
    yplot = np.maximum(yv, floor)
    low = min(xplot.min(), yplot.min())
    high = max(xplot.max(), yplot.max())
    ax.scatter(
        xplot, yplot, s=4, alpha=0.22, color=ESS,
        edgecolors="none", rasterized=True
    )
    ax.plot([low, high], [low, high], color="#333333", lw=0.9, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Official V3 probability")
    ax.set_ylabel("SpatialESS probability")
    ax.text(
        0.04, 0.96,
        "Active Jaccard = 1.000\n"
        f"max |difference| = {exact['max_abs_probability_difference']:.2e}\n"
        "p-value exact fraction = 1.000\n"
        "significance disagreements = 0",
        transform=ax.transAxes, va="top", ha="left", fontsize=7.2,
        bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "pad": 4},
    )
    clean_axes(ax, grid=False)

    fig.text(
        0.5, 0.005,
        "Ovarian Visium HD window; shared RAW-derived RCTD labels; one thread.",
        ha="center", fontsize=7, color="#555555"
    )
    save(fig, "Figure1_SPARKLE_Official_vs_SpatialESS")


def figure_raw_vs_sparkle() -> None:
    metrics_frame = pd.read_csv(DATA / "comparison_metrics.tsv", sep="\t")
    metrics = dict(zip(metrics_frame["metric"], metrics_frame["value"]))
    lr = pd.read_csv(DATA / "lr_aggregate_comparison.tsv", sep="\t")
    transition = pd.read_csv(DATA / "significant_transition.tsv", sep="\t")
    reference = pd.read_csv(DATA / "expression_reference_summary.csv")

    fig = plt.figure(figsize=(7.2, 5.9), facecolor="white")
    grid = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.34)

    ax = fig.add_subplot(grid[0, 0])
    raw_ref = reference.loc[
        reference["method"] == "RAW", "mean_snrna_corr"
    ].iloc[0]
    sparkle_ref = reference.loc[
        reference["method"] == "SPARKLE", "mean_snrna_corr"
    ].iloc[0]
    ax.set_xlim(0, 1.08)
    ax.set_ylim(-0.55, 1.55)
    panel_label(ax, "A")
    ax.set_title("Input-level effects of correction", pad=8)
    categories = ["Input counts retained", "scFFPE reference correlation"]
    raw_input = np.array([1.0, raw_ref])
    sparkle_input = np.array([0.511, sparkle_ref])
    ypos = np.array([1, 0])
    for yi, raw_value, sparkle_value in zip(ypos, raw_input, sparkle_input):
        ax.plot(
            [raw_value, sparkle_value], [yi, yi], color="#BBBBBB",
            lw=1.4, zorder=1
        )
    ax.scatter(raw_input, ypos, s=45, color=RAW, label="RAW", zorder=3)
    ax.scatter(sparkle_input, ypos, s=45, color=SPARKLE, label="SPARKLE", zorder=3)
    for yi, raw_value, sparkle_value in zip(ypos, raw_input, sparkle_input):
        ax.text(raw_value, yi + 0.16, f"{raw_value:.3f}", ha="center", fontsize=7)
        ax.text(
            sparkle_value, yi - 0.18, f"{sparkle_value:.3f}",
            ha="center", fontsize=7, color=SPARKLE, fontweight="bold"
        )
    ax.set_yticks(ypos, categories)
    ax.set_xlabel("Fraction or correlation")
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    ax.legend(frameon=False, loc="lower right", ncol=2, fontsize=7)
    clean_axes(ax, grid=False)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.55, alpha=0.75)

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    ax.set_title("Aggregate LR strength is preserved", pad=8)
    x = lr["raw_probability"].to_numpy()
    y = lr["sparkle_probability"].to_numpy()
    positive = np.concatenate([x[x > 0], y[y > 0]])
    floor = positive.min() * 0.65
    xplot = np.maximum(x, floor)
    yplot = np.maximum(y, floor)
    low = min(xplot.min(), yplot.min())
    high = max(xplot.max(), yplot.max())
    changed = y - x
    colors = np.where(changed >= 0, SPARKLE, LOST)
    ax.scatter(
        xplot, yplot, c=colors, s=8, alpha=0.48,
        edgecolors="none", rasterized=True
    )
    ax.plot([low, high], [low, high], color="#333333", lw=0.9, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("RAW aggregate LR probability")
    ax.set_ylabel("SPARKLE aggregate LR probability")
    ax.text(
        0.04, 0.96,
        f"Spearman rho = {metrics['lr_aggregate_spearman']:.3f}\n"
        f"Top-100 overlap = {int(metrics['top100_lr_overlap'])}/100",
        transform=ax.transAxes, va="top", fontsize=7.4,
        bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "pad": 4},
    )
    clean_axes(ax, grid=False)

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    ax.set_title("Communication records after correction", pad=8)
    categories = ["Active", "Significant"]
    raw_values = np.array(
        [metrics["raw_active_records"], metrics["raw_significant_records"]]
    )
    sparkle_values = np.array(
        [metrics["sparkle_active_records"], metrics["sparkle_significant_records"]]
    )
    positions = np.arange(2)
    width = 0.34
    bars_raw = ax.bar(
        positions - width / 2, raw_values, width,
        color=RAW, label="RAW", edgecolor="white"
    )
    bars_sparkle = ax.bar(
        positions + width / 2, sparkle_values, width,
        color=SPARKLE, label="SPARKLE", edgecolor="white"
    )
    ax.set_xticks(positions, categories)
    ax.set_ylabel("Communication records")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value/1000:.0f}k"))
    ax.set_ylim(0, raw_values.max() * 1.20)
    for bars in (bars_raw, bars_sparkle):
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() - raw_values.max() * 0.035,
                f"{int(bar.get_height()):,}",
                ha="center", va="top", fontsize=7, color="white",
                fontweight="bold"
            )
    active_change = 100 * (sparkle_values[0] / raw_values[0] - 1)
    sig_change = 100 * (sparkle_values[1] / raw_values[1] - 1)
    for xpos, change, ceiling in zip(positions, [active_change, sig_change], raw_values):
        ax.text(
            xpos, ceiling + raw_values.max() * 0.055, f"{change:.1f}%",
            ha="center", va="bottom", fontsize=8, color=SPARKLE,
            fontweight="bold"
        )
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    ax.set_title("Significant communication set transition", pad=8)
    labels = ["Shared", "Lost", "Gained"]
    values = transition.set_index("category").loc[
        ["Shared significant", "Lost after SPARKLE", "Gained after SPARKLE"],
        "records",
    ].to_numpy()
    colors = [SHARED, LOST, GAINED]
    bars = ax.barh(labels, values, color=colors, height=0.55)
    ax.invert_yaxis()
    ax.set_xlabel("Sender-receiver-LR records")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value/1000:.0f}k"))
    ax.set_xlim(0, values.max() * 1.28)
    for bar, value in zip(bars, values):
        ax.text(
            value + values.max() * 0.025,
            bar.get_y() + bar.get_height() / 2,
            f"{int(value):,}", va="center", fontsize=7.5, fontweight="bold"
        )
    ax.text(
        0.98, 0.08,
        f"Significant Jaccard = {metrics['significant_jaccard']:.3f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.0,
        color=SHARED, fontweight="bold"
    )
    clean_axes(ax, grid=False)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.55, alpha=0.75)

    fig.text(
        0.5, 0.005,
        "SPARKLE changes the input signal; SpatialESS inference settings are unchanged.",
        ha="center", fontsize=7, color="#555555"
    )
    save(fig, "Figure2_RAW_vs_SPARKLE_SpatialESS")


def main() -> None:
    configure()
    figure_official_vs_ess()
    figure_raw_vs_sparkle()
    print(f"Wrote publication figures to {OUT}")


if __name__ == "__main__":
    main()
