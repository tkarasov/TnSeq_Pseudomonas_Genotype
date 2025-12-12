#!/usr/bin/env bash
# Reciprocal Best Hits (RBH) pipeline using DIAMOND
# Output goes to: output_data/orthology_rbh/<timestamp>/
#
# Usage:
#   ./step0_rbh_diamond_12_2025.sh proteomeA.faa proteomeB.faa [threads]
cd /Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/scripts/orthology
set -euo pipefail

# ---------- INPUT ----------
FAA_A=DC3000_protein.faa
#"$1"
FAA_B=plate25.C2.annotation.faa
#"$2"
THREADS="${3:-4}"

# ---------- RBH THRESHOLDS ----------
EVAL_CUT="1e-10"
PID_CUT="30"
COV_CUT="0.60"

# ---------- OUTPUT DIRECTORY ----------
BASE_OUT="/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/orthology/ortholog_diamond_12_2025"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="${BASE_OUT}/rbh_${TIMESTAMP}"

mkdir -p "$OUTDIR"

echo "Output directory: $OUTDIR"

# ---------- CHECK INPUTS ----------
for f in "$FAA_A" "$FAA_B"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: Input FASTA not found or empty: $f"
    exit 1
  fi
done

command -v diamond >/dev/null || { echo "diamond not found"; exit 1; }

# ---------- TEMP DIR ----------
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

DB_A="${TMP}/A.dmnd"
DB_B="${TMP}/B.dmnd"

echo "Building DIAMOND databases..."
diamond makedb --in "$FAA_A" -d "$DB_A" >/dev/null
diamond makedb --in "$FAA_B" -d "$DB_B" >/dev/null

# Fields:
# qseqid sseqid pident length qlen slen evalue bitscore
OUTFMT="6 qseqid sseqid pident length qlen slen evalue bitscore"

RAW_A2B="${TMP}/A2B.raw"
RAW_B2A="${TMP}/B2A.raw"

echo "Running A → B..."
diamond blastp -q "$FAA_A" -d "$DB_B" -p "$THREADS" --more-sensitive -f $OUTFMT -o "$RAW_A2B"

echo "Running B → A..."
diamond blastp -q "$FAA_B" -d "$DB_A" -p "$THREADS" --more-sensitive -f $OUTFMT -o "$RAW_B2A"

# ---------- FILTER FUNCTION ----------
filter_hits () {
  IN="$1"
  OUT="$2"
  awk -v E="$EVAL_CUT" -v P="$PID_CUT" -v C="$COV_CUT" '
  BEGIN { OFS=","; print "qseqid","sseqid","pident","alen","qlen","slen","evalue","bits","cov_q","cov_s" }
  {
    q=$1; s=$2; pid=$3; alen=$4; ql=$5; sl=$6; e=$7; bits=$8;
    if (ql==0 || sl==0) next;
    covq = alen/ql;
    covs = alen/sl;
    if (e <= E && pid >= P && covq >= C && covs >= C) {
      print q, s, pid, alen, ql, sl, e, bits, covq, covs;
    }
  }' "$IN" > "$OUT"
}

FILT_A2B="${OUTDIR}/A_vs_B.filtered.csv"
FILT_B2A="${OUTDIR}/B_vs_A.filtered.csv"

echo "Filtering hits..."
filter_hits "$RAW_A2B" "$FILT_A2B"
filter_hits "$RAW_B2A" "$FILT_B2A"

# ---------- BEST HIT SELECTION ----------
best_hit () {
  IN="$1"
  OUT="$2"
  # choose best bitscore; tie by lower evalue; then higher pid
  awk -F"," '
  NR==1 { header=$0; next }
  {
    key=$1
    if (!(key in best) ||
        $8 > bestbits[key] ||
        ($8 == bestbits[key] && $7 < beste[key]) ||
        ($8 == bestbits[key] && $7 == beste[key] && $3 > bestpid[key])) {
      best[key]=$0
      bestbits[key]=$8
      beste[key]=$7
      bestpid[key]=$3
    }
  }
  END {
    print header
    for (k in best) print best[k]
  }' "$IN" | sort -t"," -k1,1 > "$OUT"
}

BEST_A2B="${OUTDIR}/A_vs_B.best.csv"
BEST_B2A="${OUTDIR}/B_vs_A.best.csv"

echo "Selecting best hits..."
best_hit "$FILT_A2B" "$BEST_A2B"
best_hit "$FILT_B2A" "$BEST_B2A"

# ---------- RECIPROCAL BEST HITS ----------
RBH="${OUTDIR}/rbh.csv"

awk -F"," '
NR==FNR && NR>1 {
  q=$1; s=$2;
  map[q]=s;    # qB → sA best hit for B→A
  next
}
NR>1 {
  a=$1; b=$2;
  if (b in map && map[b] == a) {
    print $0
  }
}' "$BEST_B2A" "$BEST_A2B" > "${RBH}.tmp"

# add header
echo "A_id,B_id,pident,alen,qlen,slen,evalue,bits,cov_q,cov_s" > "$RBH"
cat "${RBH}.tmp" >> "$RBH"
rm -f "${RBH}.tmp"

# ---------- SUMMARY ----------
SUMMARY="${OUTDIR}/rbh_summary.txt"
NUM_RBH=$(($(wc -l < "$RBH") - 1))

{
  echo "Reciprocal Best Hits Summary"
  echo "============================"
  echo "Input A: $FAA_A"
  echo "Input B: $FAA_B"
  echo "E-value cutoff: $EVAL_CUT"
  echo "Percent identity cutoff: $PID_CUT"
  echo "Coverage cutoff: $COV_CUT"
  echo ""
  echo "RBH count: $NUM_RBH"
  echo "Filtered A→B hits: $(($(wc -l < "$FILT_A2B") - 1))"
  echo "Filtered B→A hits: $(($(wc -l < "$FILT_B2A") - 1))"
  echo "Output directory: $OUTDIR"
} > "$SUMMARY"

echo ""
echo "====================================="
echo " RBH PIPELINE COMPLETE"
echo " RBH FILE: $RBH"
echo " SUMMARY: $SUMMARY"
echo "====================================="
