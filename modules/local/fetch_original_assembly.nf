process FETCH_ORIGINAL_ASSEMBLY {

  tag "${sample_id}:${assembly_accession}"
  publishDir "${params.outdir}/original_assemblies", mode: 'copy'
    conda "${projectDir}/envs/vibrio-qc.yml"

  input:
  tuple val(sample_id), val(assembly_accession)

  output:
  tuple val(sample_id),
        path("${sample_id}.original.fasta"),
        emit: assembly

  script:
  """
  datasets download genome accession \
      ${assembly_accession} \
      --include genome \
      --filename ${sample_id}.zip

  unzip -q ${sample_id}.zip -d ${sample_id}_dataset

  fasta=\$(find ${sample_id}_dataset \
      -type f \
      -name "*_genomic.fna" |
      head -n 1)

  cp "\$fasta" ${sample_id}.original.fasta
  """
}
