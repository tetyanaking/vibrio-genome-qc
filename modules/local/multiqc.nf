process MULTIQC {
    tag "cohort"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda "${projectDir}/envs/vibrio-qc.yml"
    
    input:
    path('reports/*')

    output:
    path('multiqc_report.html'), emit: report
    path('multiqc_report_data'),        emit: data

    script:
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        --config ${projectDir}/assets/multiqc_config.yaml \\
        reports/
    """
}