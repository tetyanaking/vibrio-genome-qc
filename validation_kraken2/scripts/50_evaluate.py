#!/usr/bin/env python3
"""Evaluate Kraken2 outputs against ground truth.

Reads:
  - refs/accessions.tsv  (truth taxonomy for pure runs)
  - sim/manifest.tsv     (sample -> workload/category/truth)
  - sim/mock/mock_truth.tsv (per-read truth for the mock community)
  - results/conf_<C>/<sample>.kraken.out (Kraken2 --output rows)

Writes per confidence level:
  - metrics/conf_<C>_per_sample.tsv   sens/spec/precision at genus and species
  - metrics/conf_<C>_confusion.tsv    full read-level truth x predicted taxid
  - metrics/conf_<C>_abundance.tsv    (mock only) estimated vs true proportion + L1
  - metrics/summary.tsv               aggregated across confidences

Only stdlib + pandas + NCBI taxonomy dump for genus rollup.
The taxonomy dump lives at ${TAXDUMP_DIR}/nodes.dmp / names.dmp — grab it once from
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz and unpack.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

VIBRIO_GENUS_TAXID = 662  # NCBI taxid for genus Vibrio


def load_taxonomy(taxdump_dir: Path):
    """Return (parent_of[taxid], rank_of[taxid], name_of[taxid])."""
    parent, rank, name = {}, {}, {}
    with open(taxdump_dir / "nodes.dmp") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            tid, pid, r = int(parts[0]), int(parts[1]), parts[2]
            parent[tid] = pid
            rank[tid] = r
    with open(taxdump_dir / "names.dmp") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            if parts[3] != "scientific name":
                continue
            name[int(parts[0])] = parts[1]
    return parent, rank, name


def ancestor_at_rank(taxid: int, target_rank: str, parent, rank) -> int:
    """Walk up until rank matches; 0 if none (or unclassified)."""
    if taxid == 0 or taxid not in parent:
        return 0
    seen = 0
    cur = taxid
    while cur > 1 and seen < 40:
        if rank.get(cur) == target_rank:
            return cur
        cur = parent.get(cur, 1)
        seen += 1
    return 0


def parse_kraken_output(path: Path) -> pd.DataFrame:
    """Kraken2 --output is TSV: C/U, read_id, taxid, length, kmer_string."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            rows.append((f[0], f[1], int(f[2])))
    return pd.DataFrame(rows, columns=["cu", "read_id", "pred_taxid"])


def per_read_truth_for_sample(sample_id: str, workload: str, category: str,
                              true_taxid: int, mock_truth_path: Path) -> pd.DataFrame:
    """Return DataFrame(read_id, true_taxid, true_category) for a sample."""
    if workload == "mock":
        df = pd.read_csv(mock_truth_path, sep="\t")
        return df.rename(columns={"true_taxid": "true_taxid",
                                  "true_category": "true_category"})[
            ["read_id", "true_taxid", "true_category"]
        ]
    # Pure runs: every read in the sample has the same true taxid.
    # We fill it in after we know the read ids from the Kraken output.
    return None  # sentinel — caller handles


def evaluate_one(sample: dict, kraken_dir: Path, mock_truth_path: Path,
                 parent, rank) -> dict:
    kraken_out = kraken_dir / f"{sample['sample_id']}.kraken.out"
    if not kraken_out.exists():
        return {"sample_id": sample["sample_id"], "error": "no kraken output"}

    pred = parse_kraken_output(kraken_out)
    if sample["workload"] == "mock":
        truth = pd.read_csv(mock_truth_path, sep="\t")
        truth["read_id"] = truth["read_id"].str.replace(r"/[12]$", "", regex=True)
        truth = truth.drop_duplicates(subset="read_id")
        merged = pred.merge(truth, on="read_id", how="left")
    else:
        merged = pred.copy()
        merged["true_taxid"] = int(sample["taxid"])
        merged["true_category"] = sample["category"]

    # Roll predictions + truth up to genus & species.
    merged["pred_genus"] = merged["pred_taxid"].apply(
        lambda t: ancestor_at_rank(int(t), "genus", parent, rank)
    )
    merged["true_genus"] = merged["true_taxid"].apply(
        lambda t: ancestor_at_rank(int(t), "genus", parent, rank)
    )
    merged["pred_species"] = merged["pred_taxid"].apply(
        lambda t: ancestor_at_rank(int(t), "species", parent, rank)
    )
    merged["true_species"] = merged["true_taxid"].apply(
        lambda t: ancestor_at_rank(int(t), "species", parent, rank)
    )

    # Vibrio-focused sens/spec at genus level.
    truth_vibrio = merged["true_genus"] == VIBRIO_GENUS_TAXID
    pred_vibrio = merged["pred_genus"] == VIBRIO_GENUS_TAXID
    tp = int((truth_vibrio & pred_vibrio).sum())
    fn = int((truth_vibrio & ~pred_vibrio).sum())
    fp = int((~truth_vibrio & pred_vibrio).sum())
    tn = int((~truth_vibrio & ~pred_vibrio).sum())

    sens_genus = tp / (tp + fn) if (tp + fn) else float("nan")
    spec_genus = tn / (tn + fp) if (tn + fp) else float("nan")
    prec_genus = tp / (tp + fp) if (tp + fp) else float("nan")

    # Species-level accuracy = fraction of positives whose species prediction matches truth.
    pos_reads = merged[merged["true_genus"] == VIBRIO_GENUS_TAXID]
    if len(pos_reads):
        sens_species = float(
            (pos_reads["pred_species"] == pos_reads["true_species"]).mean()
        )
    else:
        sens_species = float("nan")

    return {
        "sample_id": sample["sample_id"],
        "workload": sample["workload"],
        "category": sample["category"],
        "n_reads": len(merged),
        "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sens_genus_vibrio": round(sens_genus, 4),
        "spec_genus_vibrio": round(spec_genus, 4),
        "prec_genus_vibrio": round(prec_genus, 4),
        "acc_species_within_vibrio": round(sens_species, 4),
    }, merged


def abundance_recovery(merged_mock: pd.DataFrame) -> pd.DataFrame:
    """L1 error between true and predicted species proportions on mock."""
    true_prop = merged_mock["true_species"].value_counts(normalize=True).rename("true")
    pred_prop = merged_mock["pred_species"].value_counts(normalize=True).rename("pred")
    df = pd.concat([true_prop, pred_prop], axis=1).fillna(0.0).reset_index()
    df = df.rename(columns={"index": "species_taxid"})
    df["abs_error"] = (df["true"] - df["pred"]).abs()
    return df.sort_values("true", ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="Toolkit root")
    ap.add_argument("--taxdump", required=True,
                    help="Directory with nodes.dmp / names.dmp")
    ap.add_argument("--confidence", action="append", required=True,
                    help="Confidence value(s) to evaluate; can repeat")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    metrics_dir = root / "metrics"
    metrics_dir.mkdir(exist_ok=True)

    manifest = pd.read_csv(root / "sim" / "manifest.tsv", sep="\t").to_dict("records")
    mock_truth_path = root / "sim" / "mock" / "mock_truth.tsv"

    print("Loading NCBI taxonomy dump ...")
    parent, rank, _name = load_taxonomy(Path(args.taxdump))

    summary_rows = []
    for conf in args.confidence:
        kraken_dir = root / "results" / f"conf_{conf}"
        if not kraken_dir.exists():
            print(f"WARN: no results dir for confidence={conf}")
            continue
        per_sample = []
        mock_merged = None
        for sample in manifest:
            row, merged = evaluate_one(sample, kraken_dir, mock_truth_path, parent, rank)
            per_sample.append(row)
            if sample["workload"] == "mock":
                mock_merged = merged
            summary_rows.append({**row, "confidence": conf})

        pd.DataFrame(per_sample).to_csv(
            metrics_dir / f"conf_{conf}_per_sample.tsv", sep="\t", index=False
        )
        if mock_merged is not None:
            ab = abundance_recovery(mock_merged)
            ab.to_csv(metrics_dir / f"conf_{conf}_abundance.tsv", sep="\t", index=False)
            l1 = float(ab["abs_error"].sum())
            print(f"[conf={conf}] mock L1 abundance error = {l1:.4f}")

    pd.DataFrame(summary_rows).to_csv(
        metrics_dir / "summary.tsv", sep="\t", index=False
    )
    print(f"Wrote metrics -> {metrics_dir}")


if __name__ == "__main__":
    main()
