#!/usr/bin/env python3
"""Apply a GUNC-driven cleaning policy to a coverage/length-filtered assembly.

Reads the filtered assembly FASTA plus the per-contig coverage table and
GUNC's per-contig / per-genome contamination reports, then decides what to
do with each contig based purely on GUNC signal (FCS-GX and FCS-adaptor
are no longer part of this pipeline).

Outputs:
    <out_fasta>      cleaned assembly (contigs kept, EXCLUDEs dropped)
    <out_review>     contigs whose disposition needs a human, plus a
                     per-genome REVIEW row when GUNC failed at maxCSS
    <out_decisions>  one row per input contig with the final decision
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--coverage", required=True, type=Path)
    parser.add_argument("--gunc-summary", required=True, type=Path,
                        help="GUNC per-genome maxCSS TSV (pass.GUNC + contamination_portion)")
    parser.add_argument("--gunc-contigs", required=True, type=Path,
                        help="GUNC per-contig assignments TSV (contamination_portion per contig)")
    parser.add_argument("--gunc-exclude-threshold", type=float, default=0.30,
                        help="Per-contig contamination_portion >= this → EXCLUDE (default 0.30)")
    parser.add_argument("--gunc-review-threshold", type=float, default=0.10,
                        help="Per-contig contamination_portion in [review, exclude) → REVIEW (default 0.10)")
    parser.add_argument("--out-fasta", required=True, type=Path)
    parser.add_argument("--out-review", required=True, type=Path)
    parser.add_argument("--out-decisions", required=True, type=Path)
    return parser.parse_args()


def read_coverage(path: Path) -> dict[str, float]:
    coverage: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                coverage[row[0]] = float(row[1])
            except ValueError:
                continue
    return coverage


def read_gunc_contigs(path: Path) -> dict[str, tuple[float, str]]:
    """Map contig -> (contamination_portion, assigned_taxonomy)."""
    per_contig: dict[str, tuple[float, str]] = {}
    if not path.exists():
        return per_contig
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig = (
                row.get("contig")
                or row.get("contig_id")
                or row.get("#contig")
                or row.get("seq_id")
            )
            if not contig:
                continue
            raw_portion = (
                row.get("contamination_portion")
                or row.get("contamination_portion_maxCSS")
                or row.get("contam_portion")
                or ""
            )
            try:
                portion = float(raw_portion)
            except (TypeError, ValueError):
                continue
            taxonomy = (row.get("assigned_taxonomy") or "").strip() or "unassigned"
            per_contig[contig.strip()] = (portion, taxonomy)
    return per_contig


def read_gunc_summary(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return {k: (v or "").strip() for k, v in row.items() if k}
    return {}


def main() -> int:
    args = parse_args()

    coverage = read_coverage(args.coverage)
    gunc_contigs = read_gunc_contigs(args.gunc_contigs)
    gunc_summary = read_gunc_summary(args.gunc_summary)
    genome_failed_gunc = gunc_summary.get("pass.GUNC", "").lower() == "false"

    records = list(SeqIO.parse(str(args.assembly), "fasta"))
    if not records:
        print(f"apply_fcs_policy: {args.assembly} has no sequences", file=sys.stderr)
        return 1

    decision_columns = [
        "contig", "length", "assigned_taxonomy", "coverage",
        "gunc_contamination_portion", "resolved_action", "mode", "reason",
    ]
    review_columns = decision_columns

    kept_records = []
    decision_rows: list[dict] = []
    review_rows: list[dict] = []

    for record in records:
        contig = record.id
        cov = coverage.get(contig)
        portion, taxonomy = gunc_contigs.get(contig, (None, "unassigned"))

        base_row = {
            "contig": contig,
            "length": len(record.seq),
            "assigned_taxonomy": taxonomy,
            "coverage": "NA" if cov is None else f"{cov:.4f}",
            "gunc_contamination_portion": "NA" if portion is None else f"{portion:.4f}",
        }

        if portion is None:
            row = dict(base_row)
            row.update({
                "resolved_action": "KEEP",
                "mode": "auto",
                "reason": "no GUNC signal for contig",
            })
            decision_rows.append(row)
            kept_records.append(record)
            continue

        if portion >= args.gunc_exclude_threshold:
            row = dict(base_row)
            row.update({
                "resolved_action": "EXCLUDE",
                "mode": "auto",
                "reason": (
                    f"GUNC contamination_portion={portion:.3f} "
                    f">= {args.gunc_exclude_threshold:.2f}"
                ),
            })
            decision_rows.append(row)
            continue

        if portion >= args.gunc_review_threshold:
            row = dict(base_row)
            row.update({
                "resolved_action": "REVIEW",
                "mode": "manual",
                "reason": (
                    f"GUNC contamination_portion={portion:.3f} "
                    f"in [{args.gunc_review_threshold:.2f},"
                    f"{args.gunc_exclude_threshold:.2f})"
                ),
            })
            decision_rows.append(row)
            review_rows.append(dict(row))
            kept_records.append(record)
            continue

        row = dict(base_row)
        row.update({
            "resolved_action": "KEEP",
            "mode": "auto",
            "reason": (
                f"GUNC contamination_portion={portion:.3f} "
                f"< {args.gunc_review_threshold:.2f}"
            ),
        })
        decision_rows.append(row)
        kept_records.append(record)

    if genome_failed_gunc:
        review_rows.append({
            "contig": "__genome__",
            "length": sum(len(r.seq) for r in records),
            "assigned_taxonomy": "NA",
            "coverage": "NA",
            "gunc_contamination_portion": gunc_summary.get("contamination_portion", "NA"),
            "resolved_action": "REVIEW",
            "mode": "manual",
            "reason": (
                "GUNC pass.GUNC=False at maxCSS "
                f"(contamination_portion={gunc_summary.get('contamination_portion', 'NA')}, "
                f"n_effective_surplus_clades={gunc_summary.get('n_effective_surplus_clades', 'NA')})"
            ),
        })

    SeqIO.write(kept_records, str(args.out_fasta), "fasta")

    with args.out_decisions.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=decision_columns, delimiter="\t")
        writer.writeheader()
        for row in decision_rows:
            writer.writerow(row)

    with args.out_review.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=review_columns, delimiter="\t")
        writer.writeheader()
        for row in review_rows:
            writer.writerow(row)

    print(
        f"apply_fcs_policy: {len(records)} contigs in, "
        f"{len(kept_records)} kept, {len(review_rows)} flagged for review",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
