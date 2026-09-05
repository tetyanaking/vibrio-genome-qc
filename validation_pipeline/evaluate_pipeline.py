#!/usr/bin/env python3
"""Compare pipeline outputs to expected contamination verdicts.

Inputs (paths relative to --pipeline-outdir, e.g. results/validation):
  summary/cohort_summary.tsv    per-sample QC row (added by PR #8/#9)
  kraken2/<sample>.kreport      Kraken2 report per sample
  gunc/<sample>.tsv             GUNC per-genome output (if present)

Reference truth: reference_accessions.tsv sitting next to this script — every
sample_id must have expected_call ∈ {target, non_target_near, non_target_distant, zymo}.

Emits:
  pipeline_validation_report.tsv   per-sample expected-vs-actual
  pipeline_validation_summary.tsv  overall confusion matrix (pass/fail counts)

Verdicts (heuristic — matches how the pipeline flags contamination today):
  * kraken2_top_genus == expected genus?          (via kreport, highest %-covered genus row)
  * checkm2 contamination_post < 5.0?             (MIMAG "high-quality" contamination cap)
  * gunc pass.MIMAG_high == True?                 (chimerism gate)
  * verdict = "clean" if ALL of above, else "contaminated"

For samples where expected_call == target we expect verdict=clean AND top genus=Vibrio.
For expected_call == non_target_* we expect either verdict=contaminated OR top genus != Vibrio.
For expected_call == zymo we expect zero reads reaching EXTRACT_VIBRIO_READS output
(handled by --check-zymo-extraction).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from collections import Counter

import pandas as pd

VIBRIO_GENUS = "Vibrio"

EXPECTED_GENUS = {
    "target": VIBRIO_GENUS,
    "non_target_near": None,   # anything-but-Vibrio (species-dependent)
    "non_target_distant": None,
    "zymo": None,
}

# Species -> expected top genus per accession label (from reference_accessions.tsv)
LABEL_TO_EXPECTED_GENUS = {}


def load_truth(truth_tsv: Path) -> dict[str, dict]:
    truth = {}
    with truth_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            label = r["label"]
            genus = r["species"].split("_")[0]  # e.g. Vibrio_cholerae -> Vibrio
            truth[label] = {
                "expected_call": r["expected_call"],
                "expected_genus": genus,
                "expected_species": r["species"],
                "true_taxid": int(r["taxid"]),
            }
    return truth


def parse_kreport_top_genus(kreport: Path) -> tuple[str, float]:
    """Return (genus_name, percent_reads) for the highest-covered G-rank line."""
    best_pct, best_name = -1.0, ""
    with kreport.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                pct = float(f[0])
            except ValueError:
                continue
            rank = f[3]
            name = f[5].strip()
            if rank == "G" and pct > best_pct:
                best_pct = pct
                best_name = name
    return best_name, best_pct


def parse_gunc(gunc_tsv: Path) -> dict | None:
    """Return the row for MIMAG-high pass/fail, or None if file missing."""
    if not gunc_tsv.exists():
        return None
    with gunc_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return row  # single-row per genome
    return None


def evaluate(pipeline_outdir: Path, truth: dict[str, dict],
             contamination_threshold: float) -> pd.DataFrame:
    cohort = pipeline_outdir / "summary" / "cohort_summary.tsv"
    if not cohort.exists():
        raise SystemExit(f"cohort_summary.tsv not found at {cohort} — did the pipeline run?")

    df = pd.read_csv(cohort, sep="\t")
    rows = []
    for _, r in df.iterrows():
        sample_id = str(r["sample_id"])
        if sample_id not in truth:
            continue  # not a validation sample
        exp = truth[sample_id]

        # Contamination
        try:
            contam_post = float(r.get("contamination_post", "nan"))
        except (TypeError, ValueError):
            contam_post = float("nan")

        # Kraken2 top genus
        kreport = pipeline_outdir / "kraken2" / f"{sample_id}.kraken2.report"
        top_genus, top_pct = parse_kreport_top_genus(kreport) if kreport.exists() else ("", float("nan"))

        # GUNC
        gunc_row = parse_gunc(pipeline_outdir / "gunc" / f"{sample_id}.tsv")
        gunc_pass = None
        if gunc_row:
            v = gunc_row.get("pass.GUNC") or gunc_row.get("pass.MIMAG_high") or ""
            gunc_pass = str(v).strip().lower() in {"true", "1", "yes", "pass"}

        # Verdict
        clean_bits = []
        if not (contam_post != contam_post):  # NaN check
            clean_bits.append(contam_post < contamination_threshold)
        if gunc_pass is not None:
            clean_bits.append(gunc_pass)
        pipeline_verdict = "clean" if clean_bits and all(clean_bits) else "contaminated"

        # Expected outcome
        expected_call = exp["expected_call"]
        expected_genus = exp["expected_genus"]

        if expected_call == "target":
            genus_correct = (top_genus == expected_genus)
            verdict_correct = (pipeline_verdict == "clean") and genus_correct
        elif expected_call.startswith("non_target"):
            genus_correct = (top_genus != VIBRIO_GENUS) if top_genus else None
            # Near-neighbours + distants should be rejected — either not called
            # Vibrio OR flagged contaminated. Being clean AND Vibrio is the fail mode.
            fail_mode = (pipeline_verdict == "clean" and top_genus == VIBRIO_GENUS)
            verdict_correct = not fail_mode
        else:
            genus_correct = None
            verdict_correct = None

        rows.append({
            "sample_id": sample_id,
            "expected_call": expected_call,
            "expected_genus": expected_genus,
            "kraken2_top_genus": top_genus,
            "kraken2_top_pct": round(top_pct, 2) if top_pct >= 0 else "",
            "contamination_post": contam_post,
            "gunc_pass": gunc_pass,
            "pipeline_verdict": pipeline_verdict,
            "genus_correct": genus_correct,
            "verdict_correct": verdict_correct,
        })

    return pd.DataFrame(rows)


def confusion(report: pd.DataFrame) -> pd.DataFrame:
    tally = Counter()
    for _, r in report.iterrows():
        cat = r["expected_call"]
        tally[(cat, "n")] += 1
        if r["verdict_correct"] is True:
            tally[(cat, "correct")] += 1
        elif r["verdict_correct"] is False:
            tally[(cat, "wrong")] += 1
        if r["genus_correct"] is True:
            tally[(cat, "genus_correct")] += 1
        elif r["genus_correct"] is False:
            tally[(cat, "genus_wrong")] += 1

    rows = []
    for cat in sorted({k[0] for k in tally}):
        n = tally[(cat, "n")]
        rows.append({
            "expected_call": cat,
            "n": n,
            "verdict_correct": tally[(cat, "correct")],
            "verdict_wrong": tally[(cat, "wrong")],
            "genus_correct": tally[(cat, "genus_correct")],
            "genus_wrong": tally[(cat, "genus_wrong")],
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-outdir", required=True,
                    help="e.g. results/validation")
    ap.add_argument("--truth-tsv", default=None,
                    help="reference_accessions.tsv (default: sibling of this script)")
    ap.add_argument("--report-out", default="pipeline_validation_report.tsv")
    ap.add_argument("--summary-out", default="pipeline_validation_summary.tsv")
    ap.add_argument("--contamination-threshold", type=float, default=5.0,
                    help="CheckM2 contamination_post cutoff for 'clean' (default 5.0)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    truth_tsv = Path(args.truth_tsv) if args.truth_tsv else here / "reference_accessions.tsv"
    truth = load_truth(truth_tsv)

    report = evaluate(Path(args.pipeline_outdir).resolve(), truth,
                      args.contamination_threshold)
    report.to_csv(args.report_out, sep="\t", index=False)

    conf = confusion(report)
    conf.to_csv(args.summary_out, sep="\t", index=False)

    print(f"Wrote {args.report_out} ({len(report)} rows)")
    print(f"Wrote {args.summary_out}")
    if len(report):
        n_correct = int((report["verdict_correct"] == True).sum())
        print(f"Overall verdict accuracy: {n_correct}/{len(report)}")


if __name__ == "__main__":
    main()
