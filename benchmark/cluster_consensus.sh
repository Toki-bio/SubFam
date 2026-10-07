#!/usr/bin/env bash
# Build a SubFam-style consensus (mafft + EMBOSS cons, plurality 0.36) for every
# cluster of size >= 2 produced by another tool; singletons are kept as they are.
#   cluster_consensus.sh <seqs.fasta> <members.tsv: id TAB cluster> <out.fasta>
set -euo pipefail
FA=$1; MEM=$2; OUT=$3
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
awk -F'\t' -v dir="$W" '
    NR == FNR { cl[$1] = $2; next }
    /^>/ { id = substr($1, 2); out = dir "/" cl[id] ".fa" }
    { print > out }
' "$MEM" "$FA"
: > "$OUT"
for f in "$W"/*.fa; do
    name=$(basename "$f" .fa)
    n=$(grep -c '^>' "$f")
    if [ "$n" -lt 2 ]; then
        awk -v n="$name" '/^>/ { print ">" n; next } 1' "$f" >> "$OUT"
        continue
    fi
    plur=$(awk -v n="$n" 'BEGIN { x = n * 0.36; print (x == int(x)) ? x : int(x) + 1 }')
    mafft --thread 1 --nuc --quiet "$f" |
        cons -filter -plurality "$plur" -name "$name" |
        awk '!/^>/ { gsub(/[Nn]/, "") } 1' >> "$OUT"
done
