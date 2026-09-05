#!/usr/bin/env python3
"""
Standalone report for original (GenBank) assembly QC.

Reads cohort_summary.tsv produced by the main pipeline and writes:
  - original_qc_report.tsv  — per-sample subset of the *_original columns
                              with pass/fail flags against MIMAG thresholds
  - original_qc_report.txt  — short plaintext cohort summary

Not wired into main.nf — run manually after the pipeline finishes, e.g.

    python3 bin/original_qc_report.py results/summary/cohort_summary.tsv \\
        --outdir results/summary/original_qc

MIMAG thresholds used (Bowers et al., 2017):
  high-quality      : completeness >= 90 % and contamination < 5 %
  medium-quality    : completeness >= 50 % and contamination < 10 %
  low-quality       : anything below
BUSCO C% >= 90 is reported as an additional single-copy-marker check.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

FAILED_PREFIX = "FAILED_"

REQUIRED_COLUMNS = (
    "sample_id",
    "completeness_original",
    "contamination_original",
    "busco_C_pct_original",
)

REPORT_COLUMNS = (
    "sample_id",
    "completeness_original",
    "contamination_original",
    "busco_C_pct_original",
    "mimag_tier",
    "completeness_pass_hq",
    "contamination_pass_hq",
    "busco_pass_90",
    "notes",
)


def _to_float(value: str) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.startswith(FAILED_PREFIX):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def classify(completeness: float | None, contamination: float | None) -> str:
    if completeness is None or contamination is None:
        return "unknown"
    if completeness >= 90 and contamination < 5:
        return "high"
    if completeness >= 50 and contamination < 10:
        return "medium"
    return "low"


def build_row(row: dict[str, str]) -> dict[str, str]:
    completeness = _to_float(row.get("completeness_original", ""))
    contamination = _to_float(row.get("contamination_original", ""))
    busco = _to_float(row.get("busco_C_pct_original", ""))
    tier = classify(completeness, contamination)

    notes: list[str] = []
    if completeness is None:
        notes.append("completeness missing")
    if contamination is None:
        notes.append("contamination missing")
    if busco is None:
        notes.append("busco missing")

    return {
        "sample_id": row.get("sample_id", ""),
        "completeness_original": row.get("completeness_original", ""),
        "contamination_original": row.get("contamination_original", ""),
        "busco_C_pct_original": row.get("busco_C_pct_original", ""),
        "mimag_tier": tier,
        "completeness_pass_hq": "yes" if completeness is not None and completeness >= 90 else "no",
        "contamination_pass_hq": "yes" if contamination is not None and contamination < 5 else "no",
        "busco_pass_90": "yes" if busco is not None and busco >= 90 else "no",
        "notes": "; ".join(notes),
    }


def summarize(rows: list[dict[str, str]]) -> str:
    total = len(rows)
    tiers = {"high": 0, "medium": 0, "low": 0, "unknown": 0}
    hq_pass = 0
    busco_pass = 0
    for r in rows:
        tiers[r["mimag_tier"]] += 1
        if r["completeness_pass_hq"] == "yes" and r["contamination_pass_hq"] == "yes":
            hq_pass += 1
        if r["busco_pass_90"] == "yes":
            busco_pass += 1

    lines = [
        "Original-assembly QC cohort summary",
        "===================================",
        f"samples                   : {total}",
        f"high-quality (MIMAG)      : {tiers['high']}",
        f"medium-quality (MIMAG)    : {tiers['medium']}",
        f"low-quality (MIMAG)       : {tiers['low']}",
        f"unknown (missing data)    : {tiers['unknown']}",
        f"completeness>=90 & cont<5 : {hq_pass}",
        f"BUSCO C%>=90              : {busco_pass}",
        "",
        "Thresholds: MIMAG (Bowers et al., 2017).",
        "Source columns: completeness_original, contamination_original, busco_C_pct_original.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cohort_summary", type=Path, help="Path to cohort_summary.tsv")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Directory for output files (default: current directory).",
    )
    args = parser.parse_args(argv)

    if not args.cohort_summary.is_file():
        print(f"error: cohort summary not found: {args.cohort_summary}", file=sys.stderr)
        return 2

    with args.cohort_summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
        if missing:
            print(
                "error: cohort summary is missing required columns: "
                + ", ".join(missing)
                + "\n(re-run the pipeline with -resume after PR #9 to populate original_* columns)",
                file=sys.stderr,
            )
            return 3
        rows = [build_row(row) for row in reader]

    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.outdir / "original_qc_report.tsv"
    txt_path = args.outdir / "original_qc_report.txt"

    with tsv_path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    txt_path.write_text(summarize(rows), encoding="utf-8")

    print(f"wrote {tsv_path}")
    print(f"wrote {txt_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
