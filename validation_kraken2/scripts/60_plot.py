#!/usr/bin/env python3
"""Two plots from metrics/summary.tsv:
  1. Sensitivity vs specificity across --confidence, one point per pure-genome sample.
  2. Mock-community abundance recovery bar chart (true vs predicted proportions).

Usage: 60_plot.py --root . --confidence 0.00 0.05 0.10 0.15 0.20
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_sens_spec(summary: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    pure = summary[summary["workload"] == "pure"]
    for category, group in pure.groupby("category"):
        agg = group.groupby("confidence").agg(
            sens=("sens_genus_vibrio", "mean"),
            spec=("spec_genus_vibrio", "mean"),
        ).reset_index()
        ax.plot(1 - agg["spec"], agg["sens"], marker="o", label=category)
        for _, row in agg.iterrows():
            ax.annotate(f"c={row['confidence']}",
                        (1 - row["spec"], row["sens"]),
                        fontsize=7, textcoords="offset points", xytext=(4, 4))
    ax.set_xlabel("1 - specificity (false positive rate)")
    ax.set_ylabel("sensitivity (Vibrio recall)")
    ax.set_title("Kraken2 --confidence sweep — pure genome runs")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_abundance(ab_path: Path, out_path: Path) -> None:
    if not ab_path.exists():
        return
    ab = pd.read_csv(ab_path, sep="\t").head(15)
    x = range(len(ab))
    fig, ax = plt.subplots(figsize=(9, 5))
    w = 0.4
    ax.bar([i - w / 2 for i in x], ab["true"], width=w, label="true")
    ax.bar([i + w / 2 for i in x], ab["pred"], width=w, label="predicted")
    ax.set_xticks(list(x))
    ax.set_xticklabels(ab["species_taxid"].astype(str), rotation=45, ha="right")
    ax.set_ylabel("proportion of reads")
    ax.set_title(f"Mock community — abundance recovery ({ab_path.stem})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--confidence", nargs="+", required=True)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    metrics = root / "metrics"
    plots = root / "plots"
    plots.mkdir(exist_ok=True)

    summary = pd.read_csv(metrics / "summary.tsv", sep="\t")
    plot_sens_spec(summary, plots / "sens_vs_spec.png")

    for conf in args.confidence:
        ab_tsv = metrics / f"conf_{conf}_abundance.tsv"
        plot_abundance(ab_tsv, plots / f"abundance_conf_{conf}.png")

    print(f"Plots -> {plots}")


if __name__ == "__main__":
    main()
