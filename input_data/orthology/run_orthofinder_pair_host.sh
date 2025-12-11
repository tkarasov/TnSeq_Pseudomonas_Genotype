#!/usr/bin/env bash
# run_orthofinder_pair_host.sh
# Run OrthoFinder on two proteomes (host-only) and extract strict one-to-one pairs.
#
# Usage:
#   ./run_orthofinder_pair_host.sh <protA.faa> <protB.faa> <out_parent_dir> [threads]
#
# Example:
#   ./run_orthofinder_pair_host.sh DC3000_protein_12_2025.faa plate25.C2.annotation.faa ~/orthofinder_out 4

set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <protA.faa> <protB.faa> <out_parent_dir> [threads]"
  exit 1
fi

PROT_A="$1"
PROT_B="$2"
OUT_PARENT="$3"
THREADS="${4:-4}"

# Tools
ORTHO=$(command -v orthofinder || true)
PYTHON=$(command -v python3 || command -v python || true)
DIAMOND=$(command -v diamond || true)

if [ -z "$ORTHO" ]; then
  echo "ERROR: orthofinder not found in PATH. Please install or activate the conda env containing orthofinder."
  exit 2
fi
if [ -z "$PYTHON" ]; then
  echo "ERROR: python3 not found in PATH."
  exit 2
fi
if [ -z "$DIAMOND" ]; then
  echo "ERROR: diamond not found in PATH. OrthoFinder requires diamond for sequence search."
  exit 2
fi

# Resolve absolute paths
PROT_A=$(realpath "$PROT_A")
PROT_B=$(realpath "$PROT_B")
OUT_PARENT=$(realpath "$OUT_PARENT")

if [ ! -f "$PROT_A" ]; then echo "ERROR: protA not found: $PROT_A"; exit 1; fi
if [ ! -f "$PROT_B" ]; then echo "ERROR: protB not found: $PROT_B"; exit 1; fi
mkdir -p "$OUT_PARENT"

# Prepare temporary proteome staging dir (two files only)
TMP_PROT_DIR=$(mktemp -d /tmp/orthopair_prots_XXXXXX)
trap 'rm -rf "$TMP_PROT_DIR"' EXIT

BASENAME_A=$(basename "$PROT_A")
BASENAME_B=$(basename "$PROT_B")

# Copy into staging dir (OrthoFinder uses filenames as species IDs)
cp "$PROT_A" "$TMP_PROT_DIR/$BASENAME_A"
cp "$PROT_B" "$TMP_PROT_DIR/$BASENAME_B"

echo "Staged proteomes in: $TMP_PROT_DIR"
DATE_TAG=$(date +%Y%m%d_%H%M%S)
OUT_SUBDIR="orthofinder_pair_run_${DATE_TAG}"
OUT_DIR="${OUT_PARENT}/${OUT_SUBDIR}"
mkdir -p "$OUT_DIR"

echo "Running OrthoFinder (threads=${THREADS}) on staged proteomes..."
# Use -a 1 (analysis threads) to avoid multiprocessing spawn issues on some platforms.
# -S diamond is default; we rely on diamond being on PATH.
"$ORTHO" -f "$TMP_PROT_DIR" -t "$THREADS" -a 1 -o "$OUT_DIR" 2>&1 | tee "${OUT_DIR}/orthofinder_cli.log"

# Locate latest Results_* inside OUT_DIR (OrthoFinder places Results_* under OUT_DIR)
echo "Locating Results_* under: $OUT_DIR"
RESULTS_DIR=$(find "$OUT_DIR" -maxdepth 2 -type d -name "Results_*" -print | sort | tail -n1 || true)

if [ -z "$RESULTS_DIR" ]; then
  echo "ERROR: Could not find Results_* directory under $OUT_DIR"
  echo "Check the log: ${OUT_DIR}/orthofinder_cli.log"
  exit 3
fi

echo "Using results directory: $RESULTS_DIR"

OG_CSV="${RESULTS_DIR}/Orthogroups.csv"
if [ ! -f "$OG_CSV" ]; then
  echo "ERROR: Orthogroups.csv not found at: $OG_CSV"
  echo "Check $RESULTS_DIR for available files:"
  ls -la "$RESULTS_DIR" || true
  exit 4
fi

# Parse Orthogroups.csv to extract strict one-to-one pairs for our two species
PARSE_OUT_DIR="${RESULTS_DIR}/one2one_pairs"
mkdir -p "$PARSE_OUT_DIR"

SPECIES_A="$BASENAME_A"
SPECIES_B="$BASENAME_B"

# Run embedded Python parser
"$PYTHON" - <<PYEOF
import csv, sys, os
og_csv = os.path.expanduser("$OG_CSV")
out_dir = os.path.expanduser("$PARSE_OUT_DIR")
sp_a = "$SPECIES_A"
sp_b = "$SPECIES_B"

# Read
with open(og_csv, newline='') as fh:
    reader = csv.reader(fh)
    header = next(reader)
    species_cols = [h.strip() for h in header[1:]]
    # Normalize species names as used by OrthoFinder: filenames are used as species IDs
    # We'll ensure the given filenames match header entries
    if sp_a not in species_cols:
        print("ERROR: species A not found in Orthogroups.csv header. Available:", file=sys.stderr)
        print(", ".join(species_cols), file=sys.stderr)
        sys.exit(2)
    if sp_b not in species_cols:
        print("ERROR: species B not found in Orthogroups.csv header. Available:", file=sys.stderr)
        print(", ".join(species_cols), file=sys.stderr)
        sys.exit(2)

    orthogroups = {}
    for row in reader:
        og = row[0].strip()
        orthogroups[og] = {}
        for i, sp in enumerate(species_cols, start=1):
            cell = row[i].strip() if i < len(row) else ""
            if not cell:
                genes = []
            else:
                genes = [g.strip() for g in cell.split(';') if g.strip()]
            orthogroups[og][sp] = genes

outf = os.path.join(out_dir, f"one2one_{sp_a}_vs_{sp_b}.tsv")
with open(outf, 'w', newline='') as outfh:
    w = csv.writer(outfh, delimiter='\t')
    w.writerow(["Orthogroup","SpeciesA","GeneA","SpeciesB","GeneB"])
    for og, mapping in orthogroups.items():
        g1 = mapping.get(sp_a, [])
        g2 = mapping.get(sp_b, [])
        if len(g1) == 1 and len(g2) == 1:
            w.writerow([og, sp_a, g1[0], sp_b, g2[0]])
print("Wrote strict one-to-one pairs to:", outf)
PYEOF

echo "Done. Results (incl. one2one pairs) are in: $RESULTS_DIR"
echo "Full OrthoFinder CLI log: ${OUT_DIR}/orthofinder_cli.log"

