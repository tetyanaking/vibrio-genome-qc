process CHECKM2 {
    tag "${sample_id}_${stage}"
    publishDir "${params.outdir}/checkm2", mode: 'copy'
    conda "${projectDir}/envs/checkm2.yml"

    input:
    tuple val(sample_id), path(assembly)
    val stage

    output:
    tuple val(sample_id),
          path("${sample_id}_${stage}_checkm2/${sample_id}_${stage}_quality_report.tsv"),
          emit: quality_report
    tuple val(sample_id),
          path("${sample_id}_${stage}_checkm2"),
          emit: dir

    script:
    def db_arg = params.checkm2_db ? "--database_path ${params.checkm2_db}" : ""
    """
    mkdir -p ${sample_id}_in
    cp ${assembly} ${sample_id}_in/
    checkm2 predict \\
        --input ${sample_id}_in \\
        --output-directory ${sample_id}_${stage}_checkm2 \\
        --extension fasta \\
        --threads ${task.cpus} \\
        --force \\
        ${db_arg}
    mv ${sample_id}_${stage}_checkm2/quality_report.tsv ${sample_id}_${stage}_checkm2/${sample_id}_${stage}_quality_report.tsv
    """
}