process RESOLVE_SRA {

  tag "$sample_id"
  publishDir "${params.outdir}/metadata", mode: 'copy'
    conda "${projectDir}/envs/vibrio-qc.yml"

  input:
  tuple val(sample_id), val(assembly_accession), val(biosample)

  output:
  tuple val(sample_id),
        val(assembly_accession),
        path("${sample_id}.sra.tsv"),
        emit: runs

  script:
  """
  # Retry esearch|efetch up to 3 times to survive flaky NCBI eutils TLS drops.
  # A good runinfo dump has a header line + >=1 data row; header-only or empty
  # means eutils bailed mid-response, so back off 5s and try again.
  attempts=0
  max_attempts=3
  while [ \$attempts -lt \$max_attempts ]; do
      esearch \
          -db sra \
          -query "${biosample}[BioSample]" |
      efetch \
          -format runinfo \
          > ${sample_id}.runinfo.csv || true

      if [ -s ${sample_id}.runinfo.csv ] && [ \$(wc -l < ${sample_id}.runinfo.csv) -gt 1 ]; then
          break
      fi

      attempts=\$((attempts + 1))
      if [ \$attempts -lt \$max_attempts ]; then
          sleep 5
      fi
  done

  python3 - <<'PY'
  import csv

  with open("${sample_id}.runinfo.csv", newline="") as src, \
       open("${sample_id}.sra.tsv", "w") as out:

      reader = csv.DictReader(src)
      out.write("sample_id\\tsrr\\tplatform\\tlayout\\tstrategy\\n")

      for row in reader:
          if (
              row["Platform"] == "ILLUMINA"
              and row["LibraryLayout"] == "PAIRED"
              and row["LibraryStrategy"] == "WGS"
          ):
              out.write(
                  f"${sample_id}\\t{row['Run']}\\t"
                  f"{row['Platform']}\\t{row['LibraryLayout']}\\t"
                  f"{row['LibraryStrategy']}\\n"
              )
  PY
  """
}
