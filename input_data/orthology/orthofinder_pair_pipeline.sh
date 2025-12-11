#!/usr/bin/env bash
# orthofinder_pair_pipeline.sh
# Build docker image if needed, run OrthoFinder on two proteomes, parse Orthogroups to one2one pairs.
# Usage:
#   ./orthofinder_pair_pipeline.sh /path/to/DC3000_protein_12_2025.faa /path/to/plate25.C2.annotation.faa /path/to/host_output_dir [threads]
#
set -euo pipefail

FAA1="${1:-}"
FAA2="${2:-}"
HOST_OUT_BASE="${3:-./orthofinder_pair_out}"
THREADS="${4:-4}"

if [ -z "$FAA1" ] || [ -z "$FAA2" ]; then
  echo "Usage: $0 /path/to/proteomeA.faa /path/to/proteomeB.faa /path/to/host_output_dir [threads]" >&2
  exit 1
fi
if [ ! -f "$FAA1" ]; then echo "ERROR: FAA1 not found: $FAA1" >&2; exit 1; fi
if [ ! -f "$FAA2" ]; then echo "ERROR: FAA2 not found: $FAA2" >&2; exit 1; fi

# Normalize absolute paths
FAA1="$(cd "$(dirname "$FAA1")" && pwd)/$(basename "$FAA1")"
FAA2="$(cd "$(dirname "$FAA2")" && pwd)/$(basename "$FAA2")"
HOST_OUT_BASE="$(mkdir -p "$HOST_OUT_BASE" && cd "$HOST_OUT_BASE" && pwd)"

IMAGE_NAME="orthofinder:2.4_micromamba"
DOCKERFILE_NAME="Dockerfile.orthofinder"

# Create Dockerfile if missing
if [ ! -f "$DOCKERFILE_NAME" ]; then
  cat > "$DOCKERFILE_NAME" <<'DOCK'
FROM mambaorg/micromamba:1.4.2
SHELL ["/bin/bash", "-lc"]
RUN micromamba create -n ofenv -y -c conda-forge -c bioconda \
        orthofinder=2.4.0 \
        diamond \
        python=3.10 && \
    micromamba clean --all --yes
ENV PATH="/opt/conda/envs/ofenv/bin:${PATH}"
WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
DOCK
  echo "Wrote $DOCKERFILE_NAME"
fi

# Build image if not present
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
  echo "Building Docker image $IMAGE_NAME (this may take a few minutes)..."
  docker build -t "$IMAGE_NAME" -f "$DOCKERFILE_NAME" .
else
  echo "Docker image $IMAGE_NAME already exists."
fi

# Stage proteomes into temporary folder (host)
TMP_PROT_DIR="$(mktemp -d -t orthopair_prots_XXXXXX)"
cp -p "$FAA1" "$TMP_PROT_DIR/"
cp -p "$FAA2" "$TMP_PROT_DIR/"
echo "Proteomes staged in: $TMP_PROT_DIR"

# Create unique output dir
TS="$(date +%Y%m%d_%H%M%S)"
HOST_OUT="${HOST_OUT_BASE}/orthofinder_pair_run_${TS}"
mkdir -p "$HOST_OUT"
echo "Host output dir: $HOST_OUT"

# Run OrthoFinder in container (mount proteomes read-only, output writable)
echo "Running OrthoFinder in Docker (threads=${THREADS})..."
docker run --rm \
  -v "$TMP_PROT_DIR":/data/prot:ro \
  -v "$HOST_OUT":/data/out \
  "$IMAGE_NAME" \
  "orthofinder -f /data/prot -t ${THREADS} -o /data/out" || {
    echo "ERROR: OrthoFinder inside Docker failed. Check Docker output above and the log under $HOST_OUT/orthofinder_cli.log" >&2
    echo "Leaving staged proteomes at $TMP_PROT_DIR for debugging"
    exit 1
  }

echo "OrthoFinder finished. Results should be under: $HOST_OUT"

# find Results_* dir (most recent)
RESULTS_DIR="$(find "$HOST_OUT" -maxdepth 2 -type d -name "Results_*" -print0 | xargs -0 ls -1dt 2>/dev/null | head -n1 || true)"
if [ -z "$RESULTS_DIR" ]; then
  echo "ERROR: Could not find Results_* under $HOST_OUT" >&2
  exit 1
fi
echo "Using results directory: $RESULTS_DIR"

# find Orthogroups file
OG_FILE="$(find "$RESULTS_DIR" -maxdepth 3 -type f \( -iname "Orthogroups.csv" -o -iname "Orthogroups.tsv" \) -print | head -n1 || true)"
if [ -z "$OG_FILE" ]; then
  echo "ERROR: Orthogroups.csv/tsv not found under $RESULTS_DIR" >&2
  ls -R "$RESULTS_DIR"
  exit 1
fi
echo "Found orthogroups file: $OG_FILE"

# Determine species names used in Orthogroups header and try to match to input filenames.
SPEC_A_BASENAME="$(basename "$FAA1")"
SPEC_B_BASENAME="$(basename "$FAA2")"
# also try names without extension (common in orthofinder)
SPEC_A_NOEXT="${SPEC_A_BASENAME%.*}"
SPEC_B_NOEXT="${SPEC_B_BASENAME%.*}"

# Parse Orthogroups and write outputs
OUT_PARSE_DIR="${RESULTS_DIR}/one2one_pairs"
mkdir -p "$OUT_PARSE_DIR"
ONE2ONE_OUT="${OUT_PARSE_DIR}/one2one_${SPEC_A_BASENAME}_vs_${SPEC_B_BASENAME}.tsv"
BESTLEN_OUT="${OUT_PARSE_DIR}/best_length_pairs_${SPEC_A_BASENAME}_vs_${SPEC_B_BASENAME}.tsv"

python3 - <<PY
import csv, sys, os, json

og_file = os.path.abspath("$OG_FILE")
specA_candidates = [ "$SPEC_A_BASENAME", "$SPEC_A_NOEXT" ]
specB_candidates = [ "$SPEC_B_BASENAME", "$SPEC_B_NOEXT" ]

# Read header to get species names used by OrthoFinder
with open(og_file, 'r', newline='') as fh:
    first = fh.readline()
    # detect delimiter
    delim = ',' if ',' in first and not '\t' in first else '\t'
    header = [h.strip() for h in first.strip().split(delim)]
    species_names = header[1:]

# helper to find best matching species column name
def match_species(candidates, species_names):
    for c in candidates:
        if c in species_names:
            return c
    # try case-insensitive or partial matches
    low_names = {s.lower(): s for s in species_names}
    for c in candidates:
        if c.lower() in low_names:
            return low_names[c.lower()]
    # try prefix match
    for c in candidates:
        for s in species_names:
            if s.startswith(c) or c.startswith(s):
                return s
    return None

spA = match_species(specA_candidates, species_names)
spB = match_species(specB_candidates, species_names)

if spA is None or spB is None:
    print("ERROR: Could not match input proteome filenames to species names in Orthogroups header.", file=sys.stderr)
    print("Header species:", file=sys.stderr)
    print(", ".join(species_names), file=sys.stderr)
    print("Tried candidates for A:", specA_candidates, file=sys.stderr)
    print("Tried candidates for B:", specB_candidates, file=sys.stderr)
    sys.exit(2)

print("Matched species in Orthogroups:")
print("  A ->", spA)
print("  B ->", spB)

# reload fully with csv module
with open(og_file, newline='') as fh:
    reader = csv.reader(fh, delimiter=delim)
    header = next(reader)
    species = [s.strip() for s in header[1:]]
    ia = species.index(spA) + 1
    ib = species.index(spB) + 1

    outdir = os.path.join(os.path.dirname(og_file), "one2one_pairs")
    os.makedirs(outdir, exist_ok=True)
    one2file = os.path.join(outdir, "one2one_{}_vs_{}.tsv".format(os.path.basename("$SPEC_A_BASENAME"), os.path.basename("$SPEC_B_BASENAME")))
    bestfile = os.path.join(outdir, "best_length_pairs_{}_vs_{}.tsv".format(os.path.basename("$SPEC_A_BASENAME"), os.path.basename("$SPEC_B_BASENAME")))

    # collect lengths from input FASTA files
    def read_lengths(path):
        d={}
        name=None
        seq=[]
        with open(path) as fh2:
            for line in fh2:
                line=line.rstrip('\n')
                if not line: continue
                if line.startswith('>'):
                    if name is not None:
                        d[name]=len(''.join(seq))
                    name=line[1:].split()[0]
                    seq=[]
                else:
                    seq.append(line.strip())
            if name is not None:
                d[name]=len(''.join(seq))
        return d

    lenA = read_lengths(os.path.abspath("$FAA1"))
    lenB = read_lengths(os.path.abspath("$FAA2"))

    with open(one2file, 'w', newline='') as o1, open(bestfile, 'w', newline='') as o2:
        w1 = csv.writer(o1, delimiter='\t')
        w2 = csv.writer(o2, delimiter='\t')
        w1.writerow(["Orthogroup","SpeciesA","GeneA","SpeciesB","GeneB"])
        w2.writerow(["Orthogroup","SpeciesA","GeneA","LenA","SpeciesB","GeneB","LenB"])
        for row in reader:
            ogid = row[0].strip()
            ca = row[ia].strip() if ia < len(row) else ""
            cb = row[ib].strip() if ib < len(row) else ""
            ga = [g.strip() for g in ca.split(';') if g.strip()] if ca else []
            gb = [g.strip() for g in cb.split(';') if g.strip()] if cb else []
            if len(ga) == 1 and len(gb) == 1:
                w1.writerow([ogid, spA, ga[0], spB, gb[0]])
            if len(ga) >= 1 and len(gb) >= 1:
                # pick longest by length dict; fallback to first if length missing
                best_a = max(ga, key=lambda x: int(lenA.get(x,0)))
                best_b = max(gb, key=lambda x: int(lenB.get(x,0)))
                la = int(lenA.get(best_a,0))
                lb = int(lenB.get(best_b,0))
                w2.writerow([ogid, spA, best_a, la, spB, best_b, lb])

    print("Wrote strict one-to-one file:", one2file)
    print("Wrote best-by-length file:", bestfile)

PY

echo
echo "Pipeline complete."
echo "Results directory: $RESULTS_DIR"
echo "One-to-one and best-length files are in: $RESULTS_DIR/one2one_pairs"
echo "Temporary proteome staging dir (for debugging) left at: $TMP_PROT_DIR"
echo "Remove it when happy: rm -rf \"$TMP_PROT_DIR\""

