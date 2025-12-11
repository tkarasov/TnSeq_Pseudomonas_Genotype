#!/usr/bin/env bash
# orthofinder_pair_simple.sh
# Simple single-file pipeline: run OrthoFinder on two proteomes and extract strict one-to-one pairs.
# Usage:
#   ./orthofinder_pair_simple.sh <protA.faa> <protB.faa> <out_parent_dir> [threads]
# Example:


#./orthofinder_pair_simple.sh \
# /Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/orthology/DC3000_protein_12_2025.faa \
  /Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/orthology/plate25.C2.annotation.faa \
  /Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/orthology_host_run \
  4



set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <protA.faa> <protB.faa> <out_parent_dir> [threads]"
  exit 1
fi

PROT_A="$1"
PROT_B="$2"
OUT_PARENT="$3"
THREADS="${4:-4}"
DATE_TAG=$(date +%Y%m%d_%H%M%S)
OUT_DIR_BASE="${OUT_PARENT%/}/orthofinder_pair_${DATE_TAG}"
LOG="${OUT_DIR_BASE}/run.log"

mkdir -p "$OUT_PARENT"
mkdir -p "$OUT_DIR_BASE"
echo "Log: $LOG"
exec > >(tee -a "$LOG") 2>&1

echo "Started at $(date)"
echo "Proteome A: $PROT_A"
echo "Proteome B: $PROT_B"
echo "Host output base: $OUT_PARENT"
echo "Threads: $THREADS"

# Helper: check file exists
for f in "$PROT_A" "$PROT_B"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: file not found: $f"
    exit 2
  fi
done

BASENAME_A=$(basename "$PROT_A")
BASENAME_B=$(basename "$PROT_B")

# Try host run first
run_host() {
  echo "=== Trying host OrthoFinder run ==="
  if ! command -v orthofinder >/dev/null 2>&1; then
    echo "orthofinder not found on host PATH."
    return 1
  fi
  if ! command -v diamond >/dev/null 2>&1; then
    echo "diamond not found on host PATH (required by OrthoFinder)."
    return 1
  fi
  if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 not found on host PATH."
    return 1
  fi

  STAGE=$(mktemp -d /tmp/orthopair_host_XXXXXX)
  cp "$PROT_A" "$STAGE/$BASENAME_A"
  cp "$PROT_B" "$STAGE/$BASENAME_B"
  echo "Staged proteomes at: $STAGE"

  HOST_OUT="${OUT_DIR_BASE}/host_run"
  mkdir -p "$HOST_OUT"
  # Use -a 1 to reduce multiprocessing spawn issues on macOS; keep -t as requested for diamond threads
  echo "Running: orthofinder -f $STAGE -t $THREADS -a 1 -o $HOST_OUT"
  orthofinder -f "$STAGE" -t "$THREADS" -a 1 -o "$HOST_OUT"
  echo "Host run finished."
  # find results
  RESULTS=$(find "$HOST_OUT" -maxdepth 2 -type d -name "Results_*" -print | sort | tail -n1 || true)
  if [ -z "$RESULTS" ]; then
    echo "No Results_* found under $HOST_OUT"
    return 1
  fi
  echo "Found results: $RESULTS"
  parse_one2one "$RESULTS"
  return 0
}

# Try Docker run (if host failed)
run_docker() {
  echo "=== Trying Docker OrthoFinder run ==="
  if ! command -v docker >/dev/null 2>&1; then
    echo "docker not available."
    return 1
  fi

  IMAGE="orthofinder:2.4_fixed"   # adjust if you use another image name
  # ensure host output dir exists and is absolute
  mkdir -p "$OUT_DIR_BASE"
  STAGE=$(mktemp -d /tmp/orthopair_docker_XXXXXX)
  cp "$PROT_A" "$STAGE/$BASENAME_A"
  cp "$PROT_B" "$STAGE/$BASENAME_B"
  echo "Staged proteomes at: $STAGE"

  # map: stage -> /data/prot (ro); host out base -> /data/out_parent (rw)
  # Force PATH inside container so diamond and orthofinder (and python) are found.
  # Run container and stream output directly to host log.
  docker run --rm \
    -v "$STAGE":/data/prot:ro \
    -v "$OUT_DIR_BASE":/data/out_parent \
    --entrypoint /bin/bash \
    "$IMAGE" -lc "export PATH=/usr/bin:/opt/conda/envs/ofenv/bin:\$PATH; \
      echo 'CONTAINER PATH='\$PATH; which diamond || true; diamond version || true; which orthofinder || true; python3 --version || true; \
      orthofinder -f /data/prot -t $THREADS -a 1 -o /data/out_parent/orthofinder_run" | tee -a "$LOG"

  # check for results on host
  # two possible locations: $OUT_DIR_BASE/orthofinder_run/Results_* or $OUT_DIR_BASE/orthofinder_run/Results_*
  RESULTS=$(find "$OUT_DIR_BASE" -maxdepth 3 -type d -name "Results_*" -print | sort | tail -n1 || true)
  if [ -n "$RESULTS" ]; then
    echo "Found results: $RESULTS"
    parse_one2one "$RESULTS"
    return 0
  else
    echo "No Results_* found in $OUT_DIR_BASE after Docker run."
    return 1
  fi
}

# Parser: strict one-to-one from Orthogroups.csv
parse_one2one() {
  RESULTS_DIR="$1"
  OGCSV="${RESULTS_DIR}/Orthogroups.csv"
  if [ ! -f "$OGCSV" ]; then
    echo "ERROR: Orthogroups.csv not found at $OGCSV"
    return 1
  fi
  OUTPARSE="${RESULTS_DIR}/one2one_pairs"
  mkdir -p "$OUTPARSE"
  OUTTSV="${OUTPARSE}/one2one_${BASENAME_A}_vs_${BASENAME_B}.tsv"
  echo "Parsing $OGCSV -> $OUTTSV"

  python3 - <<PYEOF
import csv, sys, os
csvf = os.path.expanduser("${OGCSV}")
outf = os.path.expanduser("${OUTTSV}")
spA = "${BASENAME_A}"
spB = "${BASENAME_B}"
with open(csvf,newline='') as fh:
    r = csv.reader(fh)
    hdr = next(r)
    species = [h.strip() for h in hdr[1:]]
    if spA not in species or spB not in species:
        print("ERROR: species names not in header", file=sys.stderr)
        print("Header species:", ",".join(species), file=sys.stderr)
        sys.exit(2)
    rows = []
    for row in r:
        og = row[0].strip()
        mapping = {}
        for i,sp in enumerate(species, start=1):
            cell = row[i].strip() if i < len(row) else ""
            genes = [g.strip() for g in cell.split(';') if g.strip()] if cell else []
            mapping[sp] = genes
        rows.append((og,mapping))
with open(outf,'w',newline='') as outfh:
    w = csv.writer(outfh, delimiter='\t')
    w.writerow(["Orthogroup","SpeciesA","GeneA","SpeciesB","GeneB"])
    for og,m in rows:
        a = m.get(spA,[])
        b = m.get(spB,[])
        if len(a)==1 and len(b)==1:
            w.writerow([og,spA,a[0],spB,b[0]])
print("Wrote one-to-one TSV:", outf)
PYEOF

  echo "Parser complete. TSV: $OUTTSV"
  return 0
}

# Main flow: host -> docker fallback
if run_host; then
  echo "Host run + parse succeeded."
  echo "Done at $(date)"
  exit 0
else
  echo "Host run failed or incomplete. Trying Docker..."
  if run_docker; then
    echo "Docker run + parse succeeded."
    echo "Done at $(date)"
    exit 0
  else
    echo "Both host and Docker runs failed. See log: $LOG"
    echo "Done at $(date)"
    exit 5
  fi
fi
