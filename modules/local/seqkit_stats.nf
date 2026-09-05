process SEQKIT_STATS {

    tag "${sample_id}_${stage}"

    publishDir "${params.outdir}/seqkit", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

    input:
    tuple val(sample_id), path(fasta)
    val stage

    output:
    tuple val(sample_id),
          path("${sample_id}_${stage}.seqkit.stats.tsv"),
          emit: tsv

    script:
    """
    seqkit stats \
        -a \
        -T \
        ${fasta} \
        > ${sample_id}_${stage}.seqkit.stats.tsv
    """
}
