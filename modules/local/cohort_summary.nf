process COHORT_SUMMARY {

    tag "cohort"

    conda "${projectDir}/envs/vibrio-qc.yml"

    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path summary_tsvs

    output:
    path("cohort_summary.tsv"), emit: tsv

    script:
    """
    python3 - <<'PY'
    from pathlib import Path

    files = sorted(Path(".").glob("*.summary.tsv"))
    if not files:
        raise SystemExit("cohort_summary: no per-sample summary TSVs found")

    header_written = False
    with open("cohort_summary.tsv", "w", encoding="utf-8") as out:
        for path in files:
            with path.open(encoding="utf-8") as handle:
                lines = handle.read().splitlines()

            if not lines:
                continue

            if not header_written:
                out.write(lines[0] + "\\n")
                header_written = True

            for line in lines[1:]:
                if line.strip():
                    out.write(line + "\\n")
    PY
    """
}
