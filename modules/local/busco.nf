process BUSCO {

    tag "${sample_id}"
    publishDir "${params.outdir}/busco", mode: 'copy'
    conda "${projectDir}/envs/busco.yml"

    input:
    tuple val(sample_id), path(assembly)
    val lineage
    val stage_tag

    output:
    tuple val(sample_id),
          path("${sample_id}_${stage_tag}_busco/short_summary*.txt"),
          emit: summary

    tuple val(sample_id),
          path("${sample_id}_${stage_tag}_busco"),
          emit: dir

    script:
    """
    busco \
        --in ${assembly} \
        --lineage_dataset ${lineage} \
        --mode genome \
        --out ${sample_id}_${stage_tag}_busco \
        --cpu ${task.cpus} \
        --offline \
        --download_path \${BUSCO_DOWNLOAD_PATH}
    """
}
