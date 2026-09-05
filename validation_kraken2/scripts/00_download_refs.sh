#!/usr/bin/env bash
# Download RefSeq reference genomes listed in refs/accessions.tsv.
# Requires: NCBI `datasets` CLI (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
# Output: refs/fasta/<label>.fna
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TSV="${ROOT}/refs/accessions.tsv"
OUT="${ROOT}/refs/fasta"
mkdir -p "${OUT}"

command -v datasets >/dev/null 2>&1 || {
    echo "ERROR: NCBI 'datasets' CLI not found on PATH." >&2
    echo "Install: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/" >&2
    exit 1
}

tail -n +2 "${TSV}" | while IFS=$'\t' read -r acc category species taxid label; do
    dest="${OUT}/${label}.fna"
    if [[ -s "${dest}" ]]; then
        echo "[skip] ${label} — already downloaded"
        continue
    fi
    echo "[get ] ${acc} -> ${dest}"
    tmp="$(mktemp -d)"
    datasets download genome accession "${acc}" \
        --include genome \
        --filename "${tmp}/pkg.zip" >/dev/null
    unzip -q "${tmp}/pkg.zip" -d "${tmp}/unpack"
    fna="$(find "${tmp}/unpack" -name '*_genomic.fna' -o -name '*.fna' | head -n1)"
    if [[ -z "${fna}" ]]; then
        echo "ERROR: no FASTA in package for ${acc}" >&2
        exit 1
    fi
    # Rewrite headers so every contig carries an unambiguous label — makes downstream
    # ground-truth mapping trivial (read id will inherit the header prefix).
    awk -v L="${label}" -v T="${taxid}" '
        /^>/ { n++; print ">" L "_ctg" n " taxid=" T; next }
        { print }
    ' "${fna}" > "${dest}"
    rm -rf "${tmp}"
done

echo "Done. FASTAs in ${OUT}"
