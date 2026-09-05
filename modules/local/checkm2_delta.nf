process CHECKM2_DELTA {

    tag "${sample_id}"

    conda "${projectDir}/envs/vibrio-qc.yml"

    publishDir "${params.outdir}/checkm2_delta", mode: 'copy'

    input:
    tuple val(sample_id), path(pre_report, stageAs: 'pre_quality_report.tsv'), path(post_report, stageAs: 'post_quality_report.tsv')

    output:
    tuple val(sample_id),
          path("${sample_id}.checkm2_delta.tsv"),
          emit: tsv

    script:
    """
    python3 - <<'PY'
    import csv
    from pathlib import Path

    sample_id = "${sample_id}"
    pre_path = Path("pre_quality_report.tsv")
    post_path = Path("post_quality_report.tsv")
    out_path = Path("${sample_id}.checkm2_delta.tsv")


    def first_row(path):
        with path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            for row in reader:
                return row
        return {}


    def to_float(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None


    pre = first_row(pre_path)
    post = first_row(post_path)

    completeness_pre = to_float(pre.get("Completeness"))
    completeness_post = to_float(post.get("Completeness"))
    contamination_pre = to_float(pre.get("Contamination"))
    contamination_post = to_float(post.get("Contamination"))


    def delta(before, after):
        if before is None or after is None:
            return ""
        return f"{after - before:.4f}"


    def formatted(value):
        if value is None:
            return ""
        return f"{value:.4f}"


    columns = [
        "sample_id",
        "completeness_pre",
        "completeness_post",
        "completeness_delta",
        "contamination_pre",
        "contamination_post",
        "contamination_delta",
    ]

    row = [
        sample_id,
        formatted(completeness_pre),
        formatted(completeness_post),
        delta(completeness_pre, completeness_post),
        formatted(contamination_pre),
        formatted(contamination_post),
        delta(contamination_pre, contamination_post),
    ]

    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("\\t".join(columns) + "\\n")
        handle.write("\\t".join(row) + "\\n")
    PY
    """
}
