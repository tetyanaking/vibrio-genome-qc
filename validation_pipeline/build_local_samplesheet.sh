#!/usr/bin/env bash
set -euo pipefail

OUT="local_samplesheet.csv"
echo "sample_id,assembly_accession,reads_r1,reads_r2" > "$OUT"

tail -n +2 data/accessions.tsv | while IFS=$'\t' read -r accession category species taxid label; do
    r1=$(find sim/pure -iname "${label}_R1.fastq.gz" | head -1)
    r2=$(find sim/pure -iname "${label}_R2.fastq.gz" | head -1)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: missing reads for ${label} (accession ${accession}) — skipping" >&2
        continue
    fi

    echo "${label},${accession},$(realpath "$r1"),$(realpath "$r2")" >> "$OUT"
done

echo
echo "Wrote ${OUT}:"
cat "$OUT"
