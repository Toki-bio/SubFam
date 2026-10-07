#!/usr/bin/env bash
# Partition with vsearch at a loose identity, run SubFam (k-mer order) inside every partition,
# pool partitions smaller than N into one "rest" partition. Evaluate like the benchmark.
#   partition_test.sh <prefix> <id> [N] [threads]     (prefix.fasta, prefix.truth.tsv, prefix.masters.fasta)
set -euo pipefail
P=$1; ID=$2; N=${3:-50}; T=${4:-2}; B=$(cd "$(dirname "$P")" && pwd); S=$(basename "$P")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
vsearch --cluster_fast "$B/$S.fasta" --id "$ID" --threads "$T" --quiet --uc "$W/p.uc"
awk -F'\t' '$1=="S"{print $9"\t"$2} $1=="H"{print $9"\t"$2}' "$W/p.uc" > "$W/part.tsv"
# cluster sizes; small ones -> "rest"
awk -F'\t' -v n="$N" 'NR==FNR{c[$2]++; next} {print $1"\t"(c[$2]>=n ? "p"$2 : "rest")}' "$W/part.tsv" "$W/part.tsv" > "$W/part2.tsv"
cut -f2 "$W/part2.tsv" | sort | uniq -c | sort -k1,1nr | awk '{printf "%s:%s ", $2, $1} END{print ""}' | cut -c1-150
mkdir -p "$W/parts"
awk -F'\t' -v d="$W/parts" 'NR==FNR{p[$1]=$2; next} /^>/{id=substr($1,2); f=d"/"p[id]".fasta"} {print > f}' "$W/part2.tsv" "$B/$S.fasta"
: > "$W/all.cons.fasta"; : > "$W/all.chunks.tsv"
for f in "$W"/parts/*.fasta; do
    p=$(basename "$f" .fasta)
    "$(dirname "$0")/../SubFam.sh" -n "$N" -t "$T" -o "$W/sf_$p" -x "$p" "$f" > /dev/null
    cat "$W/sf_$p/$p.cons.fasta" >> "$W/all.cons.fasta"
    cat "$W/sf_$p/$p.chunks.tsv" >> "$W/all.chunks.tsv"
done
python3 "$(dirname "$0")/evaluate.py" "vsearch $ID + SubFam -n $N" "$W/all.cons.fasta" "$W/all.chunks.tsv" "$B/$S.truth.tsv" "$B/$S.masters.fasta"
