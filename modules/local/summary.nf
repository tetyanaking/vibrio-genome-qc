process SAMPLE_SUMMARY {

    tag "${sample_id}"

    conda "${projectDir}/envs/vibrio-qc.yml"

    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    tuple val(sample_id),
          path(seqkit_stats_original, stageAs: 'original_seqkit_stats.tsv'),
          path(checkm2_original),
          path(busco_original, stageAs: 'original_busco_summary.txt'),
          path(seqkit_stats, stageAs: 'post_seqkit_stats.tsv'),
          path(checkm2_delta),
          path(busco_summary, stageAs: 'post_busco_summary.txt'),
          path(fcs_decisions)

    output:
    tuple val(sample_id),
          path("${sample_id}.summary.tsv"),
          emit: tsv

    script:
    """
    python3 - <<'PY'
    import csv
    import re
    from pathlib import Path


    sample_id = "${sample_id}"
    NA_MARKER = "NA.placeholder"
    FAILED = "FAILED_"


    def missing(path):
        return Path(path).name == NA_MARKER


    def first_row(path):
        with Path(path).open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            for row in reader:
                return row
        return {}


    def busco_c_pct(path, is_missing):
        if is_missing:
            return FAILED
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        match = re.search(r"C:\\s*([0-9.]+)%", text)
        return match.group(1) if match else ""


    stats_original_missing   = missing("original_seqkit_stats.tsv")
    checkm2_original_missing = missing("${checkm2_original}")
    busco_original_missing   = missing("original_busco_summary.txt")
    stats_missing = missing("post_seqkit_stats.tsv")
    delta_missing = missing("${checkm2_delta}")
    busco_missing = missing("post_busco_summary.txt")
    fcs_missing   = missing("${fcs_decisions}")

    stats_original   = {} if stats_original_missing else first_row("original_seqkit_stats.tsv")
    checkm2_original = {} if checkm2_original_missing else first_row("${checkm2_original}")
    stats = {} if stats_missing else first_row("post_seqkit_stats.tsv")
    delta = {} if delta_missing else first_row("${checkm2_delta}")

    busco_c_original = busco_c_pct("original_busco_summary.txt", busco_original_missing)
    busco_c = busco_c_pct("post_busco_summary.txt", busco_missing)

    if fcs_missing:
        exclude_count = trim_count = review_count = FAILED
    else:
        exclude_count = trim_count = review_count = 0
        with Path("${fcs_decisions}").open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            for row in reader:
                action = (row.get("resolved_action") or row.get("decision") or "").strip().upper()
                if action == "EXCLUDE":
                    exclude_count += 1
                elif action == "TRIM":
                    trim_count += 1
                elif action == "REVIEW":
                    review_count += 1


    def s(source_missing, value):
        return FAILED if source_missing else value


    columns = [
        "sample_id",
        "completeness_original",
        "contamination_original",
        "busco_C_pct_original",
        "n_contigs",
        "total_len",
        "N50",
        "completeness_pre",
        "completeness_post",
        "completeness_delta",
        "contamination_pre",
        "contamination_post",
        "contamination_delta",
        "busco_C_pct",
        "fcs_exclude",
        "fcs_trim",
        "fcs_review",
    ]

    row = [
        sample_id,
        s(checkm2_original_missing, checkm2_original.get("Completeness", "")),
        s(checkm2_original_missing, checkm2_original.get("Contamination", "")),
        busco_c_original,
        s(stats_missing, stats.get("num_seqs", "")),
        s(stats_missing, stats.get("sum_len", "")),
        s(stats_missing, stats.get("N50", "")),
        s(delta_missing, delta.get("completeness_pre", "")),
        s(delta_missing, delta.get("completeness_post", "")),
        s(delta_missing, delta.get("completeness_delta", "")),
        s(delta_missing, delta.get("contamination_pre", "")),
        s(delta_missing, delta.get("contamination_post", "")),
        s(delta_missing, delta.get("contamination_delta", "")),
        busco_c,
        str(exclude_count),
        str(trim_count),
        str(review_count),
    ]

    with open(f"{sample_id}.summary.tsv", "w", encoding="utf-8") as out:
        out.write("\\t".join(columns) + "\\n")
        out.write("\\t".join(row) + "\\n")
    PY
    """
}
