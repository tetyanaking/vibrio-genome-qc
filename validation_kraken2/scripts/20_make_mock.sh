#!/usr/bin/env bash
# Build a mock community by sub-sampling the pure-genome FASTQs at fixed proportions.
# Composition (default): 60% Vibrio (positives) / 20% near / 20% distant.
# Output: sim/mock/mock_R{1,2}.fastq.gz + sim/mock/mock_truth.tsv (per-read taxid).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TSV="${ROOT}/refs/accessions.tsv"
PURE="${ROOT}/sim/pure"
OUT="${ROOT}/sim/mock"
MANIFEST="${ROOT}/sim/manifest.tsv"
mkdir -p "${OUT}"

TOTAL_PAIRS="${MOCK_PAIRS:-1500000}"
POS_FRAC="${MOCK_POS:-0.60}"
NEAR_FRAC="${MOCK_NEAR:-0.20}"
DIST_FRAC="${MOCK_DIST:-0.20}"

command -v seqtk >/dev/null 2>&1 || {
    echo "ERROR: seqtk not found (needed for sub-sampling). Install with conda: conda install -c bioconda seqtk" >&2
    exit 1
}

# Count members per category
declare -A CAT_COUNT
while IFS=$'\t' read -r acc category species taxid label; do
    [[ "${acc}" == "accession" ]] && continue
    CAT_COUNT[${category}]=$(( ${CAT_COUNT[${category}]:-0} + 1 ))
done < "${TSV}"

R1_OUT="${OUT}/mock_R1.fastq.gz"
R2_OUT="${OUT}/mock_R2.fastq.gz"
TRUTH="${OUT}/mock_truth.tsv"

: > "${OUT}/.r1.list"
: > "${OUT}/.r2.list"
printf "read_id\ttrue_taxid\ttrue_label\ttrue_category\n" > "${TRUTH}"

tail -n +2 "${TSV}" | while IFS=$'\t' read -r acc category species taxid label; do
    r1="${PURE}/${label}_R1.fastq.gz"
    r2="${PURE}/${label}_R2.fastq.gz"
    [[ -s "${r1}" ]] || { echo "WARN: missing ${r1}"; continue; }

    case "${category}" in
        positive) frac="${POS_FRAC}";  members="${CAT_COUNT[positive]}"  ;;
        near)     frac="${NEAR_FRAC}"; members="${CAT_COUNT[near]}"      ;;
        distant)  frac="${DIST_FRAC}"; members="${CAT_COUNT[distant]}"   ;;
        *) continue ;;
    esac

    n=$(python3 -c "import math; print(int(math.floor(${TOTAL_PAIRS} * ${frac} / ${members})))")
    echo "[mock] ${label} — ${n} pairs"

    tmp_r1="${OUT}/.sub_${label}_R1.fastq"
    tmp_r2="${OUT}/.sub_${label}_R2.fastq"
    # seed = hashed label -> reproducible + independent per genome
    seed=$(python3 -c "import zlib,sys; print(zlib.crc32(sys.argv[1].encode()) % 100000)" "${label}")
    seqtk sample -s "${seed}" "${r1}" "${n}" > "${tmp_r1}"
    seqtk sample -s "${seed}" "${r2}" "${n}" > "${tmp_r2}"

    # Truth rows — one per read pair; read id is the first token after '@'
    awk -v L="${label}" -v T="${taxid}" -v C="${category}" '
        NR%4==1 { sub(/^@/,""); split($0,a," "); print a[1] "\t" T "\t" L "\t" C }
    ' "${tmp_r1}" >> "${TRUTH}"

    cat "${tmp_r1}" >> "${OUT}/.mock_R1.fastq"
    cat "${tmp_r2}" >> "${OUT}/.mock_R2.fastq"
    rm -f "${tmp_r1}" "${tmp_r2}"
done

echo "[gzip] compressing mock FASTQs"
gzip -c "${OUT}/.mock_R1.fastq" > "${R1_OUT}"
gzip -c "${OUT}/.mock_R2.fastq" > "${R2_OUT}"
rm -f "${OUT}/.mock_R1.fastq" "${OUT}/.mock_R2.fastq" "${OUT}/.r1.list" "${OUT}/.r2.list"

# Register mock in manifest (append)
if ! grep -q $'\tmock\t' "${MANIFEST}"; then
    printf "mock_community\tmock\tmixed\tmixed\t0\tmixed\t%s\t%s\t%d\n" \
        "${R1_OUT}" "${R2_OUT}" "${TOTAL_PAIRS}" >> "${MANIFEST}"
fi

echo "Mock ready: ${R1_OUT}"
echo "Truth: ${TRUTH}"
