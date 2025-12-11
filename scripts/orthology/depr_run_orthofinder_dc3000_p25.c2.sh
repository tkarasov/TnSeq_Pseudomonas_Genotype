#!/usr/bin/env bash
# run_orthofinder_dc3000_p25.c2.sh
# Forces OrthoFinder CLI to run single-threaded (works around macOS/conda multiprocessing.spawn errors).
# Usage:
#   ./run_orthofinder_force_singlethread.sh <faa1> <faa2> [out_base]
# Defaults use your repo paths.


set -euo pipefail

FAA1="${1:-/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/orthology/DC3000_protein_12_2025.faa}"
FAA2="${2:-/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/orthology/plate25.C2.annotation.faa}"
OUT_BASE="${3:-/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_pair_cli}"

# checks
if [ ! -f "$FAA1" ]; then echo "ERROR: FAA1 not found: $FAA1" >&2; exit 1; fi
if [ ! -f "$FAA2" ]; then echo "ERROR: FAA2 not found: $FAA2" >&2; exit 1; fi

# find orthofinder CLI
ORTH_CMD="$(command -v orthofinder || true)"
if [ -z "$ORTH_CMD" ]; then
  ORTH_CMD="/Users/talia/opt/anaconda3/bin/orthofinder"
  if [ ! -x "$ORTH_CMD" ]; then
    echo "ERROR: orthofinder not found in PATH and fallback not executable: $ORTH_CMD" >&2
    exit 1
  fi
fi
echo "Using OrthoFinder CLI: $ORTH_CMD"

# stage two proteomes in a temporary dir
TMP_PROT_DIR="$(mktemp -d -t orthofinder_pair_XXXXXX)"
trap 'rm -rf "$TMP_PROT_DIR"' EXIT
cp -p "$FAA1" "$TMP_PROT_DIR/"
cp -p "$FAA2" "$TMP_PROT_DIR/"
echo "Staged proteomes in: $TMP_PROT_DIR"

SPEC_A="$(basename "$FAA1")"
SPEC_B="$(basename "$FAA2")"

# unique output dir so OrthoFinder won't refuse to run
TS="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_BASE}_run_${TS}"
mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/orthofinder_cli.log"

echo "Running OrthoFinder single-threaded (this avoids multiprocessing spawn errors)..."
# macOS helpful env vars: OBJC_ var can help some fork-related issues
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
# ensure OrthoFinder uses 1 thread
if ! "$ORTH_CMD" -f "$TMP_PROT_DIR" -t 1 -o "$OUT_DIR" >"$LOG" 2>&1 ; then
  echo "ERROR: OrthoFinder CLI failed. Showing log head ($LOG):" >&2
  sed -n '1,200p' "$LOG" >&2 || true
  exit 1
fi

echo "OrthoFinder completed. Log: $LOG"
echo "Results should be in: $OUT_DIR (look for Results_* directory)"
# print the Results_* candidate(s)
ls -ld "$OUT_DIR"/Results_* 2>/dev/null || true
echo "Done."