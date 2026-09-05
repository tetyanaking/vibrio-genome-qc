#!/usr/bin/env bash
# Classify every FASTQ pair in the manifest against k2_standard_16_GB_20260626.
# Output: results/conf_<C>/<sample_id>.kraken.out and .kreport
# Usage: 30_run_kraken2.sh [confidence]        (default 0.10)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${ROOT}/sim/manifest.tsv"

: "${KRAKEN2_DB:?Set KRAKEN2_DB to the path of your k2_standard_16_GB_20260626 directory}"
CONF="${1:-0.10}"
THREADS="${K2_THREADS:-8}"
OUT="${ROOT}/results/conf_${CONF}"
mkdir -p "${OUT}"

command -v kraken2 >/dev/null 2>&1 || {
    echo "ERROR: kraken2 not on PATH" >&2; exit 1
}

echo "DB=${KRAKEN2_DB}   confidence=${CONF}   threads=${THREADS}"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample workload label species taxid category r1 r2 nreads; do
    out="${OUT}/${sample}.kraken.out"
    rep="${OUT}/${sample}.kreport"
    if [[ -s "${out}" ]]; then
        echo "[skip] ${sample}"
        continue
    fi
    echo "[k2  ] ${sample} (${workload}, ${category})"
    kraken2 \
        --db "${KRAKEN2_DB}" \
        --threads "${THREADS}" \
        --confidence "${CONF}" \
        --paired \
        --gzip-compressed \
        --output "${out}" \
        --report "${rep}" \
        "${r1}" "${r2}" >/dev/null
done

echo "Done. Outputs in ${OUT}"
