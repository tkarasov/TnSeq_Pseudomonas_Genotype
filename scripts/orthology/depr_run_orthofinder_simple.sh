#!/usr/bin/env bash
# Usage:
#   ./run_orthofinder_simple.sh <proteomeA.faa> <proteomeB.faa> <output_directory>

set -euo pipefail

PROT_A="$1"
PROT_B="$2"
OUT_PARENT="$3"

# --- Create output directory -----------------------------------------
mkdir -p "$OUT_PARENT"
DATE_TAG=$(date +%Y%m%d_%H%M%S)
OUT_RUN="${OUT_PARENT}/orthofinder_run_${DATE_TAG}"
mkdir -p "$OUT_RUN"

echo "Running OrthoFinder in Docker..."
echo "Proteome A: $PROT_A"
echo "Proteome B: $PROT_B"
echo "Output dir: $OUT_RUN"

# --- Stage the two proteomes for OrthoFinder --------------------------
STAGE_DIR=$(mktemp -d /tmp/orthopair_XXXXXX)
cp "$PROT_A" "$STAGE_DIR/$(basename "$PROT_A")"
cp "$PROT_B" "$STAGE_DIR/$(basename "$PROT_B")"

# --- Run OrthoFinder inside Docker ------------------------------------
docker run --rm \
  -v "$STAGE_DIR":/data/prot:ro \
  -v "$OUT_RUN":/data/out_parent \
  orthofinder:2.4_fixed \
  /bin/bash -lc "export PATH=/usr/bin:/opt/conda/envs/ofenv/bin:\$PATH; \
                 which diamond; diamond version; \
                 which orthofinder; \
                 orthofinder -f /data/prot -t 4 -o /data/out_parent/Results"

echo "OrthoFinder run complete."

# --- Locate Orthogroups.tsv -------------------------------------------
OG_FILE="$(find ${OUT_RUN}/Results -maxdepth 3 -type f -name 'Orthogroups.tsv' -print | head -n1)"

if [ -z "$OG_FILE" ]; then
  echo "ERROR: Could not find Orthogroups.tsv in OrthoFinder output."
  exit 1
fi

echo "Found Orthogroups file:"
echo "  $OG_FILE"

# --- Parse one-to-one orthologs ---------------------------------------
OUTDIR="${OUT_RUN}/one2one"
mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR}/one2one_$(basename "$PROT_A")_vs_$(basename "$PROT_B").tsv"

echo "Extracting strict one-to-one orthologs..."
python3 - <<PY
import csv, sys, os

og = "$OG_FILE"
protA = os.path.basename("$PROT_A")
protB = os.path.basename("$PROT_B")
out = "$OUTFILE"

with open(og, newline='') as fh:
    reader = csv.reader(fh, delimiter='\t')
    header = next(reader)
    species = header[1:]
    # match exact or extension-stripped names
    def matches(p, s):
        p0 = p.split('.')[0]
        s0 = s.split('.')[0]
        return p == s or p0 == s0 or p0 == s.replace('_','-') or p0 == s.replace('_','')

    tryA = [s for s in species if matches(protA, s)]
    tryB = [s for s in species if matches(protB, s)]
    if not tryA or not tryB:
        print("Species names in Orthogroups.tsv do not match proteome basenames.")
        print("Header species:", species)
        sys.exit(1)

    spA = tryA[0]
    spB = tryB[0]

    ia = species.index(spA) + 1
    ib = species.index(spB) + 1

    with open(out, 'w', newline='') as outfh:
        w = csv.writer(outfh, delimiter='\t')
        w.writerow(["Orthogroup","SpeciesA","GeneA","SpeciesB","GeneB"])
        for row in reader:
            ogid = row[0]
            genesA = [g for g in row[ia].split(';') if g.strip()]
            genesB = [g for g in row[ib].split(';') if g.strip()]
            if len(genesA) == 1 and len(genesB) == 1:
                w.writerow([ogid, spA, genesA[0], spB, genesB[0]])

print("Wrote one-to-one TSV:")
print(out)
PY

echo ""
echo "Done."
echo "Output directory:"
echo "  $OUT_RUN"
echo "One-to-one ortholog file:"
echo "  $OUTFILE"

