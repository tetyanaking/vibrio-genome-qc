process MERGE_RUNS {

  tag "$sample_id"
  publishDir "${params.outdir}/merged_reads", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id),
        path("${sample_id}_R{1,2}.fastq.gz"),
        emit: reads

  script:
  def r1 = reads.findAll { it.name.endsWith('_1.fastq.gz') }
  def r2 = reads.findAll { it.name.endsWith('_2.fastq.gz') }

  if (!r1 || !r2) {
        error "Missing paired FASTQ files for ${sample_id}: R1=${r1.size()}, R2=${r2.size()}"
  }

  """
  cat ${r1.join(' ')} > ${sample_id}_R1.fastq.gz
  cat ${r2.join(' ')} > ${sample_id}_R2.fastq.gz
  """
}
