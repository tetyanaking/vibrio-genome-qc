process FASTP {

    tag "${sample_id}"

    label 'medium'

    publishDir "${params.outdir}/fastp", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_R{1,2}.trimmed.fastq.gz"),
          emit: reads

    tuple val(sample_id),
          path("${sample_id}.fastp.json"),
          emit: json

    tuple val(sample_id),
          path("${sample_id}.fastp.html"),
          emit: html

    script:
    """
    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --length_required 50 \
        --json ${sample_id}.fastp.json \
        --html ${sample_id}.fastp.html \
        --thread ${task.cpus}
    """
}
