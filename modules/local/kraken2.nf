process KRAKEN2 {

    tag "$sample_id"

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

    input:
    tuple val(sample_id), path(reads)
    path kraken2_db

    output:
    tuple val(sample_id),
          path("${sample_id}.kraken2.out"),
          emit: classified

    tuple val(sample_id),
          path("${sample_id}.kraken2.report"),
          emit: report

    script:
    """
    kraken2 \
        --db ${kraken2_db} \
        --paired ${reads[0]} ${reads[1]} \
        --output ${sample_id}.kraken2.out \
        --report ${sample_id}.kraken2.report \
        --threads ${task.cpus} \
        --memory-mapping \
	    --use-names
    """
}
