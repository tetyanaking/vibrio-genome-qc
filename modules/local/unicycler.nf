process UNICYCLER {
  tag "$sample_id"
  publishDir "${params.outdir}/unicycler", mode: 'copy'
    conda "${projectDir}/envs/vibrio-qc.yml"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.scaffolds.fasta"), emit: scaffolds
  tuple val(sample_id), path("${sample_id}.assembly.gfa"),    emit: gfa
  tuple val(sample_id), path("${sample_id}.unicycler.log"),   emit: log

  script:
  """
  unicycler \\
      -1 ${reads[0]} -2 ${reads[1]} \\
      -o ${sample_id}_unicycler \\
      -t ${task.cpus} \\
      --mode normal

  cp ${sample_id}_unicycler/assembly.fasta ${sample_id}.scaffolds.fasta
  cp ${sample_id}_unicycler/assembly.gfa   ${sample_id}.assembly.gfa
  cp ${sample_id}_unicycler/unicycler.log  ${sample_id}.unicycler.log
  """
}
