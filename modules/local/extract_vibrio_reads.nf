process EXTRACT_VIBRIO_READS {

  tag "$sample_id"

  publishDir "${params.outdir}/vibrio_reads", mode: 'copy'

    conda "${projectDir}/envs/vibrio-qc.yml"

  input:
  tuple val(sample_id), path(reads), path(kraken_out), path(kraken_report)
  val vibrio_taxid

  output:
  tuple val(sample_id),
        path("${sample_id}_vibrio_R{1,2}.fastq.gz"),
        emit: reads

  tuple val(sample_id),
        path("${sample_id}.extract.stats.tsv"),
        emit: stats

  script:
  """
  extract_kraken_reads.py \
      -k ${kraken_out} \
      -r ${kraken_report} \
      -s1 ${reads[0]} \
      -s2 ${reads[1]} \
      -o ${sample_id}_vibrio_R1.fastq \
      -o2 ${sample_id}_vibrio_R2.fastq \
      -t ${vibrio_taxid} \
      --include-children \
      --fastq-output

  pigz -p ${task.cpus} \
      ${sample_id}_vibrio_R1.fastq \
      ${sample_id}_vibrio_R2.fastq

  n_in=\$(zcat ${reads[0]} | awk 'NR%4==1' | wc -l)
  n_out=\$(zcat ${sample_id}_vibrio_R1.fastq.gz | awk 'NR%4==1' | wc -l)

  pct=\$(awk -v a=\$n_out -v b=\$n_in \
      'BEGIN {if (b==0) printf "0.00"; else printf "%.2f", 100*a/b}')

  printf "sample_id\\tread_pairs_in\\tread_pairs_vibrio\\tpct\\n" \
      > ${sample_id}.extract.stats.tsv

  printf "%s\\t%d\\t%d\\t%s\\n" \
      "${sample_id}" "\$n_in" "\$n_out" "\$pct" \
      >> ${sample_id}.extract.stats.tsv
  """
}