process GUNC {

    tag "${sample_id}"

    conda "${projectDir}/envs/gunc.yml"

    publishDir "${params.outdir}/gunc", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path gunc_db

    output:
    tuple val(sample_id),
          path("${sample_id}.gunc.tsv"),
          emit: summary

    tuple val(sample_id),
          path("${sample_id}.gunc.contig_assignments.tsv"),
          emit: contigs

    tuple val(sample_id),
          path("diamond_output/*.diamond.out"),
          emit: diamond, optional: true

    script:
    """
    mkdir -p input
    cp ${assembly} input/${sample_id}.fa

    gunc run \\
        --input_dir input \\
        --db_file ${gunc_db} \\
        --out_dir . \\
        --file_suffix .fa \\
        --detailed_output \\
        --threads ${task.cpus}

    # Normalise the maxCSS-level summary filename
    mv GUNC.*.maxCSS_level.tsv ${sample_id}.gunc.tsv

    # Normalise the per-contig assignments filename. Always write the
    # target file so downstream joins never break; write a header-only
    # placeholder if GUNC did not emit per-contig data (small, clean genome).
    contig_asg=\$(find gene_counts -name '*.contig_assignments.tsv' -print -quit 2>/dev/null || true)
    if [ -n "\${contig_asg}" ] && [ -f "\${contig_asg}" ]; then
        mv "\${contig_asg}" ${sample_id}.gunc.contig_assignments.tsv
    else
        printf "contig\\tassigned_taxonomy\\tcontamination_portion\\n" \\
            > ${sample_id}.gunc.contig_assignments.tsv
    fi
    """

    stub:
    """
    printf "genome\\tn_contigs_flagged\\tcontamination_portion\\tpass.GUNC\\n" \\
        > ${sample_id}.gunc.tsv
    printf "%s\\t0\\t0.00\\tTrue\\n" ${sample_id} \\
        >> ${sample_id}.gunc.tsv

    printf "contig\\tassigned_taxonomy\\tcontamination_portion\\n" \\
        > ${sample_id}.gunc.contig_assignments.tsv
    """
}
