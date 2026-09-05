#!/usr/bin/env bash
# Simulate paired-end Illumina reads from every reference in refs/fasta.
# Requires: InSilicoSeq (`pip install InSilicoSeq`).
# Output: sim/pure/<label>_R{1,2}.fastq.gz plus sim/manifest.tsv (ground truth).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTA_DIR="${ROOT}/refs/fasta"
TSV="${ROOT}/refs/accessions.tsv"
OUT="${ROOT}/sim/pure"
MANIFEST="${ROOT}/sim/manifest.tsv"
mkdir -p "${OUT}"

MODEL="${ISS_MODEL:-miseq}"        # miseq | hiseq | novaseq | or a trained .npz
N_READS="${ISS_READS:-200000}"     # per-genome pairs
CPUS="${ISS_CPUS:-4}"

command -v iss >/dev/null 2>&1 || {
    echo "ERROR: InSilicoSeq (iss) not found. Install: pip install InSilicoSeq" >&2
    exit 1
}

printf "sample_id\tworkload\tsource_label\tspecies\ttaxid\tcategory\tfastq_r1\tfastq_r2\ttrue_reads\n" > "${MANIFEST}"

tail -n +2 "${TSV}" | while IFS=$'\t' read -r acc category species taxid label; do
    fna="${FASTA_DIR}/${label}.fna"
    prefix="${OUT}/${label}"
    r1="${prefix}_R1.fastq.gz"
    r2="${prefix}_R2.fastq.gz"

    if [[ ! -s "${fna}" ]]; then
        echo "WARN: missing FASTA for ${label}, skipping" >&2
        continue
    fi

    if [[ -s "${r1}" && -s "${r2}" ]]; then
        echo "[skip] ${label} — reads already simulated"
    else
        echo "[sim ] ${label} — ${N_READS} pairs, model=${MODEL}"
        iss generate \
            --genomes "${fna}" \
            --model "${MODEL}" \
            --n_reads "${N_READS}" \
            --cpus "${CPUS}" \
            --output "${prefix}" \
            --compress \
            --quiet
        # iss writes _R1.fastq.gz / _R2.fastq.gz already
    fi

    printf "%s\tpure\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n" \
        "pure_${label}" "${label}" "${species}" "${taxid}" "${category}" \
        "${r1}" "${r2}" "${N_READS}" >> "${MANIFEST}"
done

echo "Manifest: ${MANIFEST}"
