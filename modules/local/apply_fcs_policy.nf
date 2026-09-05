process APPLY_FCS_POLICY {

    tag "${sample_id}"

    publishDir "${params.outdir}/fcs_cleaned", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

    input:
    tuple val(sample_id),
          path(filtered_fa),
          path(cov_tsv),
          path(gunc_summary),
          path(gunc_contigs)

    output:
    tuple val(sample_id),
          path("${sample_id}.clean.fasta"),
          emit: clean_assembly

    tuple val(sample_id),
          path("${sample_id}.review.tsv"),
          emit: review

    tuple val(sample_id),
          path("${sample_id}.decisions.tsv"),
          emit: decisions

    script:
    """
    python3 "${projectDir}/bin/apply_fcs_policy.py" \\
        --assembly "${filtered_fa}" \\
        --coverage "${cov_tsv}" \\
        --gunc-summary "${gunc_summary}" \\
        --gunc-contigs "${gunc_contigs}" \\
        --out-fasta "${sample_id}.clean.fasta" \\
        --out-review "${sample_id}.review.tsv" \\
        --out-decisions "${sample_id}.decisions.tsv"
    """

    stub:
    """
    cp "${filtered_fa}" "${sample_id}.clean.fasta"

    printf "contig\\tlength\\tassigned_taxonomy\\tcoverage\\tgunc_contamination_portion\\tresolved_action\\tmode\\treason\\n" \\
        > "${sample_id}.review.tsv"

    printf "contig\\tlength\\tassigned_taxonomy\\tcoverage\\tgunc_contamination_portion\\tresolved_action\\tmode\\treason\\n" \\
        > "${sample_id}.decisions.tsv"
    """
}
