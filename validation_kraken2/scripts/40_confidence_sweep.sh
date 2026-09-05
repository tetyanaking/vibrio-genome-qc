#!/usr/bin/env bash
# Re-run classification at multiple --confidence thresholds so you can pick an
# operating point from a sensitivity-vs-specificity curve instead of guessing.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SWEEP="${K2_SWEEP:-0.00 0.05 0.10 0.15 0.20}"
for c in ${SWEEP}; do
    echo "=== confidence=${c} ==="
    bash "${ROOT}/scripts/30_run_kraken2.sh" "${c}"
done
