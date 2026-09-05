process COV_LENGTH_FILTER {

    tag "${sample_id}"

    label 'medium'

    publishDir "${params.outdir}/coverage_filter", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

    input:
    tuple val(sample_id),
          path(assembly),
          path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.filtered.fasta"),
          emit: filtered

    tuple val(sample_id),
          path("${sample_id}.cov.tsv"),
          emit: coverage

    tuple val(sample_id),
          path("${sample_id}.removed.fasta"),
          emit: removed

    tuple val(sample_id),
          path("${sample_id}.contig_decisions.tsv"),
          emit: decisions

    tuple val(sample_id),
          path("${sample_id}.mapped.bam"),
          path("${sample_id}.mapped.bam.bai"),
          emit: bam

    script:
    if (reads.size() != 2) {
        error """
        COV_LENGTH_FILTER expected exactly two paired-end read files for
        ${sample_id}, but received: ${reads}
        """.stripIndent()
    }

    def r1 = reads[0]
    def r2 = reads[1]

    """
    set -euo pipefail

    minimap2 \
        -ax sr \
        -t ${task.cpus} \
        "${assembly}" \
        "${r1}" \
        "${r2}" |
    samtools sort \
        -@ ${task.cpus} \
        -o "${sample_id}.mapped.bam" \
        -

    samtools index "${sample_id}.mapped.bam"

    samtools depth \
        -aa \
        "${sample_id}.mapped.bam" |
    awk '
        {
            sum[\$1] += \$3
            count[\$1] += 1
        }
        END {
            print "contig\\tmean_coverage"
            for (contig in sum) {
                print contig "\\t" sum[contig] / count[contig]
            }
        }
    ' > "${sample_id}.cov.tsv"

    python3 - <<'PY'
    from pathlib import Path
    import sys

    import numpy as np
    from Bio import SeqIO


    sample_id = "${sample_id}"
    assembly_path = Path("${assembly}")
    coverage_path = Path("${sample_id}.cov.tsv")

    filtered_path = Path("${sample_id}.filtered.fasta")
    removed_path = Path("${sample_id}.removed.fasta")
    decisions_path = Path("${sample_id}.contig_decisions.tsv")

    minimum_length = 500
    lower_coverage_factor = 0.25
    upper_coverage_factor = 3.0


    # Read mean coverage table.
    coverage = {}

    with coverage_path.open(encoding="utf-8") as handle:
        header = next(handle, None)

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\\n").split("\\t")

            if len(fields) != 2:
                raise ValueError(
                    f"Malformed coverage row at line {line_number}: {line!r}"
                )

            contig_id, mean_coverage = fields

            try:
                coverage[contig_id] = float(mean_coverage)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid coverage value at line {line_number}: "
                    f"{mean_coverage!r}"
                ) from exc


    positive_coverage = [
        value
        for value in coverage.values()
        if np.isfinite(value) and value > 0
    ]

    if not positive_coverage:
        raise RuntimeError(
            f"No positive contig coverage values were found for {sample_id}"
        )


    # Median is being used as a robust estimate of the main coverage peak.
    # It is not a true statistical mode.
    central_coverage = float(np.median(positive_coverage))

    minimum_coverage = lower_coverage_factor * central_coverage
    maximum_coverage = upper_coverage_factor * central_coverage


    kept_records = []
    removed_records = []
    decision_rows = []

    assembly_records = list(
        SeqIO.parse(str(assembly_path), "fasta")
    )

    if not assembly_records:
        raise RuntimeError(
            f"No sequences were found in assembly: {assembly_path}"
        )


    for record in assembly_records:
        contig_id = record.id
        length = len(record.seq)
        mean_coverage = coverage.get(contig_id)

        reasons = []

        if length < minimum_length:
            reasons.append("length_below_500")

        if mean_coverage is None:
            reasons.append("coverage_missing")
        else:
            if mean_coverage < minimum_coverage:
                reasons.append("coverage_below_0.25x_median")

            if mean_coverage > maximum_coverage:
                reasons.append("coverage_above_3x_median")

        if reasons:
            decision = "remove"
            removed_records.append(record)
        else:
            decision = "keep"
            kept_records.append(record)

        decision_rows.append(
            {
                "contig": contig_id,
                "length": length,
                "mean_coverage": (
                    "NA"
                    if mean_coverage is None
                    else f"{mean_coverage:.4f}"
                ),
                "central_coverage": f"{central_coverage:.4f}",
                "minimum_coverage": f"{minimum_coverage:.4f}",
                "maximum_coverage": f"{maximum_coverage:.4f}",
                "decision": decision,
                "reason": (
                    ";".join(reasons)
                    if reasons
                    else "passed_length_and_coverage_filters"
                ),
            }
        )


    if not kept_records:
        raise RuntimeError(
            f"Filtering removed every contig for {sample_id}. "
            "Review the coverage thresholds before continuing."
        )


    SeqIO.write(
        kept_records,
        str(filtered_path),
        "fasta",
    )

    SeqIO.write(
        removed_records,
        str(removed_path),
        "fasta",
    )


    columns = [
        "contig",
        "length",
        "mean_coverage",
        "central_coverage",
        "minimum_coverage",
        "maximum_coverage",
        "decision",
        "reason",
    ]

    with decisions_path.open("w", encoding="utf-8") as handle:
        handle.write("\\t".join(columns) + "\\n")

        for row in decision_rows:
            handle.write(
                "\\t".join(str(row[column]) for column in columns)
                + "\\n"
            )


    print(
        f"{sample_id}: central coverage = {central_coverage:.2f}x; "
        f"retained {len(kept_records)} of {len(assembly_records)} contigs",
        file=sys.stderr,
    )
    PY
    """

    stub:
    """
    printf ">${sample_id}_contig_1\\nACGTACGTACGT\\n" \
        > "${sample_id}.filtered.fasta"

    printf ">${sample_id}_removed_contig\\nACGT\\n" \
        > "${sample_id}.removed.fasta"

    printf "contig\\tmean_coverage\\n${sample_id}_contig_1\\t50.0\\n" \
        > "${sample_id}.cov.tsv"

    printf "contig\\tlength\\tmean_coverage\\tcentral_coverage\\tminimum_coverage\\tmaximum_coverage\\tdecision\\treason\\n" \
        > "${sample_id}.contig_decisions.tsv"

    touch "${sample_id}.mapped.bam"
    touch "${sample_id}.mapped.bam.bai"
    """
}
