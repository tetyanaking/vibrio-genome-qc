process DOWNLOAD_READS {

  tag "${sample_id}:${srr}"
  publishDir "${params.outdir}/raw_reads", mode: 'copy'
    conda "${projectDir}/envs/vibrio-qc.yml"

  cpus 4

  input:
  tuple val(sample_id), val(srr)

  output:
  tuple val(sample_id),
        val(srr),
        path("${sample_id}.${srr}_{1,2}.fastq.gz"),
        emit: reads

  script:
  """
  prefetch ${srr} --max-size u

  fasterq-dump \
      ${srr}/${srr}.sra \
      --split-files \
      --threads ${task.cpus}

  mv ${srr}_1.fastq ${sample_id}.${srr}_1.fastq
  mv ${srr}_2.fastq ${sample_id}.${srr}_2.fastq

  pigz \
      -p ${task.cpus} \
      ${sample_id}.${srr}_1.fastq \
      ${sample_id}.${srr}_2.fastq
  """
}
