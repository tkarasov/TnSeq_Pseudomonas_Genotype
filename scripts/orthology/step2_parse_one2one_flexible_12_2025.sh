#!/usr/bin/env bash
# parse_one2one_flexible.sh
# Finds the Orthogroups file (tsv) and extracts strict one-to-one pairs for two proteomes,
# trying many heuristics to match the proteome names to the Orthogroups header.
#
# Usage:
#   ./parse_one2one_flexible.sh <Orthogroups.tsv> <protA.faa> <protB.faa>
# Example:
#   ./parse_one2one_flexible.sh /path/to/Orthogroups.tsv DC3000_protein_12_2025.faa plate25.C2.annotation.faa

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <Orthogroups.tsv> <protA.faa> <protB.faa>"
  exit 1
fi

OG_TSV="$1"
PROT_A="$2"
PROT_B="$3"

if [ ! -f "$OG_TSV" ]; then
  echo "ERROR: Orthogroups file not found: $OG_TSV"
  exit 2
fi

OUTDIR="$(dirname "$OG_TSV")/one2one_pairs"
mkdir -p "$OUTDIR"

python3 - <<PY
import csv, os, sys, re
og = os.path.expanduser("$OG_TSV")
protA = os.path.basename("$PROT_A")
protB = os.path.basename("$PROT_B")
def cand_names(p):
    # produce many candidate species names from a filename
    names = set()
    names.add(p)
    p_noext = re.sub(r'(\.faa|\.fa|\.fasta|\.pep|\.faa.gz|\.fa.gz)$','',p, flags=re.I)
    names.add(p_noext)
    names.add(p_noext.replace('.','_'))
    names.add(p_noext.replace('.',''))
    names.add(p_noext.replace('-','_'))
    names.add(p_noext.replace('.','_').replace('-','_'))
    # also try only the first token before first dot or underscore
    names.add(re.split(r'[._-]', p_noext)[0])
    # try uppercase/lowercase variants
    names.add(p_noext.lower())
    names.add(p_noext.upper())
    return [n for n in names if n]

def find_best(header, candidates):
    # header: list of species names (trimmed)
    # candidates: list of possible names for one proteome
    header_norm = [h.strip() for h in header]
    # exact match (case-sensitive)
    for c in candidates:
        if c in header_norm:
            return c
    # exact match case-insensitive
    for c in candidates:
        for h in header_norm:
            if c.lower() == h.lower():
                return h
    # header contains candidate or candidate contains header (longer/shorter)
    for c in candidates:
        for h in header_norm:
            if c in h or h in c:
                return h
    # replace dots/underscores and compare
    def norm(x): return re.sub(r'[\._\-]','',x).lower()
    for c in candidates:
        for h in header_norm:
            if norm(c) == norm(h):
                return h
    # startswith
    for c in candidates:
        for h in header_norm:
            if h.lower().startswith(c.lower()) or c.lower().startswith(h.lower()):
                return h
    # no match
    return None

with open(og, newline='') as fh:
    reader = csv.reader(fh, delimiter='\t')
    try:
        hdr = next(reader)
    except StopIteration:
        print("ERROR: empty Orthogroups file", file=sys.stderr); sys.exit(3)
    species = [h.strip() for h in hdr[1:]]
    # produce candidate names
    candA = cand_names(protA)
    candB = cand_names(protB)
    matchA = find_best(species, candA)
    matchB = find_best(species, candB)
    if not matchA or not matchB:
        print("ERROR: Could not map proteome filenames to Orthogroups header.", file=sys.stderr)
        print("Proteome A candidates (first 10):", candA[:10], file=sys.stderr)
        print("Proteome B candidates (first 10):", candB[:10], file=sys.stderr)
        print("Header species (full):", file=sys.stderr)
        print(", ".join(species), file=sys.stderr)
        sys.exit(4)
    # If both resolve to the same header name (ambiguous), warn and fail
    if matchA == matchB:
        print("ERROR: Both proteomes resolved to same species name in header:", matchA, file=sys.stderr)
        print("Header species:", ", ".join(species), file=sys.stderr)
        sys.exit(5)
    # We will re-run parsing now using the detected header names
    outfn = os.path.join(os.path.dirname(og), "one2one_pairs", f"one2one_{protA}_vs_{protB}.tsv")
    # rewind file and parse full
    fh.seek(0)
    reader = csv.reader(fh, delimiter='\t')
    hdr = next(reader)
    species = [h.strip() for h in hdr[1:]]
    rows = []
    for row in reader:
        ogid = row[0].strip()
        mapping = {}
        for i,sp in enumerate(species, start=1):
            cell = row[i].strip() if i < len(row) else ""
            genes = [g.strip() for g in cell.split(';') if g.strip()] if cell else []
            mapping[sp] = genes
        rows.append((ogid,mapping))
    # write TSV
    import csv
    with open(outfn, 'w', newline='') as outfh:
        w = csv.writer(outfh, delimiter='\t')
        w.writerow(["Orthogroup","SpeciesA","GeneA","SpeciesB","GeneB"])
        hits = 0
        for ogid,m in rows:
            a = m.get(matchA,[])
            b = m.get(matchB,[])
            if len(a)==1 and len(b)==1:
                w.writerow([ogid, matchA, a[0], matchB, b[0]])
                hits += 1
    print("Wrote TSV:", outfn)
    print("Proteome A matched header name:", matchA)
    print("Proteome B matched header name:", matchB)
    print("Number of strict one-to-one pairs written:", hits)
PY

