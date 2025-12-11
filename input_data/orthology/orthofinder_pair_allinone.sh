#!/usr/bin/env bash
# orthofinder_pair_allinone.sh (fixed)
# Build image if needed, run OrthoFinder on two proteomes, parse orthogroups -> TSVs.
# Usage:
#   ./orthofinder_pair_allinone.sh /path/to/proteomeA.faa /path/to/proteomeB.faa /path/to/host_output_parent [threads]
#
set -euo pipefail
IFS=$'\n\t'

FAA1="${1:-}"
FAA2="${2:-}"
HOST_OUT_PARENT="${3:-./orthofinder_pair_out}"
THREADS="${4:-4}"

if [ -z "$FAA1" ] || [ -z "$FAA2" ]; then
  echo "Usage: $0 /path/to/proteomeA.faa /path/to/proteomeB.faa /path/to/host_output_parent [threads]" >&2
  exit 1
fi
if [ ! -f "$FAA1" ]; then echo "ERROR: FAA1 not found: $FAA1" >&2; exit 1; fi
if [ ! -f "$FAA2" ]; then echo "ERROR: FAA2 not found: $FAA2" >&2; exit 1; fi

# normalize absolute paths
FAA1="$(cd "$(dirname "$FAA1")" && pwd)/$(basename "$FAA1")"
FAA2="$(cd "$(dirname "$FAA2")" && pwd)/$(basename "$FAA2")"
HOST_OUT_PARENT="$(mkdir -p "$HOST_OUT_PARENT" && cd "$HOST_OUT_PARENT" && pwd)"

IMAGE_NAME="orthofinder:2.4_fixed"
DOCKERFILE_NAME="Dockerfile.orthofinder.clean"

echo "FAA1: $FAA1"
echo "FAA2: $FAA2"
echo "Host output parent: $HOST_OUT_PARENT"
echo "Threads: $THREADS"
echo "Docker image: $IMAGE_NAME"

# Create Dockerfile (idempotent)
cat > "$DOCKERFILE_NAME" <<'DOCK'
# Dockerfile.orthofinder.clean
FROM mambaorg/micromamba:1.4.2
SHELL ["/bin/bash", "-lc"]
RUN micromamba create -n ofenv -y -c conda-forge -c bioconda \
        orthofinder=2.4.0 \
        diamond \
        python=3.10 && \
    micromamba clean --all --yes
ENV CONDA_PREFIX=/opt/conda/envs/ofenv
ENV PATH="${CONDA_PREFIX}/bin:${PATH}"
RUN if [ -x "${CONDA_PREFIX}/bin/python" ]; then \
      ln -sf "${CONDA_PREFIX}/bin/python" /usr/bin/python3 || true; \
    fi
WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
DOCK

# Build image if not present
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
  echo "Building Docker image $IMAGE_NAME (may take a few minutes)..."
  docker build -t "$IMAGE_NAME" -f "$DOCKERFILE_NAME" .
else
  echo "Docker image $IMAGE_NAME already exists, skipping build."
fi

# Stage proteomes into a temporary directory on host
TMP_PROT_DIR="$(mktemp -d -t orthopair_prots_XXXXXX)"
cp -p "$FAA1" "$TMP_PROT_DIR/"
cp -p "$FAA2" "$TMP_PROT_DIR/"
echo "Proteomes copied to temporary staging dir: $TMP_PROT_DIR"

# Create a unique subdirectory name, but DO NOT create it on the host — let OrthoFinder create it
TS="$(date +%Y%m%d_%H%M%S)"
OUT_SUBDIR="orthofinder_pair_run_${TS}"
echo "OrthoFinder will create subdir: $OUT_SUBDIR under $HOST_OUT_PARENT"

# Run OrthoFinder inside container using env python (explicit path) and pass a non-existent -o subdir
echo "Running OrthoFinder in Docker (streaming logs to host)..."
docker run --rm \
  -v "$TMP_PROT_DIR":/data/prot:ro \
  -v "$HOST_OUT_PARENT":/data/out_parent \
  "$IMAGE_NAME" \
  "/opt/conda/envs/ofenv/bin/python /opt/conda/envs/ofenv/bin/orthofinder -f /data/prot -t ${THREADS} -o /data/out_parent/${OUT_SUBDIR}" \
  2>&1 | tee "${HOST_OUT_PARENT}/${OUT_SUBDIR}_cli_partial.log" || {
    echo "ERROR: OrthoFinder inside Docker failed. Check log: ${HOST_OUT_PARENT}/${OUT_SUBDIR}_cli_partial.log" >&2
    echo "Staged proteomes left at: $TMP_PROT_DIR" >&2
    exit 1
  }

# After successful run, the output directory should now exist under host parent
HOST_OUT="${HOST_OUT_PARENT}/${OUT_SUBDIR}"
echo "OrthoFinder finished. Host output dir: $HOST_OUT"
echo "Full run log: $HOST_OUT/orthofinder_cli.log (if produced) and the streamed partial log above."

# locate Results_* directory (most recent) under the created HOST_OUT
RESULTS_DIR="$(find "$HOST_OUT" -maxdepth 2 -type d -name "Results_*" -print0 | xargs -0 ls -1dt 2>/dev/null | head -n1 || true)"
if [ -z "$RESULTS_DIR" ]; then
  echo "ERROR: Could not find Results_* under $HOST_OUT" >&2
  ls -R "$HOST_OUT"
  exit 1
fi
echo "Using results directory: $RESULTS_DIR"

# find Orthogroups.csv (or tsv)
OG_FILE="$(find "$RESULTS_DIR" -maxdepth 3 -type f \( -iname "Orthogroups.csv" -o -iname "Orthogroups.tsv" \) -print | head -n1 || true)"
if [ -z "$OG_FILE" ]; then
  echo "ERROR: Orthogroups.csv/tsv not found under $RESULTS_DIR" >&2
  ls -R "$RESULTS_DIR"
  exit 1
fi
echo "Found orthogroups file: $OG_FILE"

# Derive species identifiers from input filenames
SPEC_A_BASENAME="$(basename "$FAA1")"
SPEC_B_BASENAME="$(basename "$FAA2")"
SPEC_A_NOEXT="${SPEC_A_BASENAME%.*}"
SPEC_B_NOEXT="${SPEC_B_BASENAME%.*}"

# Parse Orthogroups -> produce one2one and best-length TSVs (run on host python)
OUT_PARSE_DIR="${RESULTS_DIR}/one2one_pairs"
mkdir -p "$OUT_PARSE_DIR"

ONE2ONE_OUT="${OUT_PARSE_DIR}/one2one_${SPEC_A_BASENAME}_vs_${SPEC_B_BASENAME}.tsv"
BESTLEN_OUT="${OUT_PARSE_DIR}/best_length_pairs_${SPEC_A_BASENAME}_vs_${SPEC_B_BASENAME}.tsv"

python3 - <<PY
import csv, os, sys

og_file = os.path.abspath("${OG_FILE}")
specA_candidates = [ "${SPEC_A_BASENAME}", "${SPEC_A_NOEXT}" ]
specB_candidates = [ "${SPEC_B_BASENAME}", "${SPEC_B_NOEXT}" ]

with open(og_file, 'r', newline='') as fh:
    first = fh.readline()
    if '\t' in first and ',' not in first:
        delim = '\t'
    elif ',' in first and '\t' not in first:
        delim = ','
    else:
        delim = ',' if ',' in first else '\t'
    header = [h.strip() for h in first.strip().split(delim)]
    species_names = header[1:]

def match_species(candidates, species_names):
    for c in candidates:
        if c in species_names:
            return c
    lower = {s.lower(): s for s in species_names}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    for c in candidates:
        for s in species_names:
            if s.startswith(c) or c.startswith(s) or c in s or s in c:
                return s
    return None

spA = match_species(specA_candidates, species_names)
spB = match_species(specB_candidates, species_names)

if spA is None or spB is None:
    print("ERROR: Could not match input filenames to species in Orthogroups header.", file=sys.stderr)
    print("Header species:", file=sys.stderr)
    print(", ".join(species_names), file=sys.stderr)
    print("Candidates A:", specA_candidates, file=sys.stderr)
    print("Candidates B:", specB_candidates, file=sys.stderr)
    sys.exit(2)

print("Matched species names:")
print("  A ->", spA)
print("  B ->", spB)

with open(og_file, newline='') as fh:
    reader = csv.reader(fh, delimiter=delim)
    header = next(reader)
    species = [s.strip() for s in header[1:]]
    ia = species.index(spA) + 1
    ib = species.index(spB) + 1

    def read_lengths(path):
        d = {}
        name = None
        seq_parts = []
        with open(path) as fh2:
            for line in fh2:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if name is not None:
                        d[name] = len(''.join(seq_parts))
                    name = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if name is not None:
                d[name] = len(''.join(seq_parts))
        return d

    lenA = read_lengths(os.path.abspath("${FAA1}"))
    lenB = read_lengths(os.path.abspath("${FAA2}"))

    one2file = os.path.join("${OUT_PARSE_DIR}", "one2one_{}_vs_{}.tsv".format(os.path.basename("${SPEC_A_BASENAME}"), os.path.basename("${SPEC_B_BASENAME}")))
    bestfile = os.path.join("${OUT_PARSE_DIR}", "best_length_pairs_{}_vs_{}.tsv".format(os.path.basename("${SPEC_A_BASENAME}"), os.path.basename("${SPEC_B_BASENAME}")))

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
                best_a = max(ga, key=lambda x: int(lenA.get(x,0)))
                best_b = max(gb, key=lambda x: int(lenB.get(x,0)))
                la = int(lenA.get(best_a,0))
                lb = int(lenB.get(best_b,0))
                w2.writerow([ogid, spA, best_a, la, spB, best_b, lb])

    print("Wrote strict one-to-one file:", one2file)
    print("Wrote best-by-length file:", bestfile)

PY

echo
echo "Pipeline finished."
echo "Results (orthogroups) directory: $RESULTS_DIR"
echo "Parsed outputs are in: $OUT_PARSE_DIR"
echo "Strict one-to-one file: $ONE2ONE_OUT"
echo "Best-by-length file: $BESTLEN_OUT"
echo
echo "Temporary proteome staging dir (left for debugging): $TMP_PROT_DIR"
echo "Delete it when happy: rm -rf \"$TMP_PROT_DIR\""

