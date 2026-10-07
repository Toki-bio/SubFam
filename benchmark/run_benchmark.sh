#!/usr/bin/env bash
# SubFam vs. identity-threshold clustering (CD-HIT, VSEARCH, MMseqs2) on simulated
# repeat families with known subfamilies.
#   run_benchmark.sh [outdir] [threads]
# Requires: mafft, EMBOSS, cd-hit-est, vsearch, mmseqs, python3.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-bench_out}; T=${2:-$(nproc)}
mkdir -p "$OUT"; cd "$OUT"

now() { date +%s.%N; }
secs() { awk -v a="$1" -v b="$(now)" 'BEGIN { printf "%.1f", b - a }'; }

run() {  # run <label> <reps.fasta> <members.tsv> <start time>
    printf '%s\t%s\n' "$(python3 "$HERE/evaluate.py" "$1" "$2" "$3" "$S.truth.tsv" "$S.masters.fasta")" "$(secs "$4")"
}

printf 'scenario\tmethod\treps\treps_ge10\tsingleton_frac\tpurity\tmedian_id\tsubfam_recovered\tseconds\n' > results.tsv

for scen in young:0.03 middle:0.08 old:0.15; do
    S=${scen%%:*}
    python3 "$HERE/simulate.py" "$S" --divergence "${scen#*:}"
    {
    # SubFam
    for n in 50 20; do
        t=$(now)
        "$HERE/../SubFam.sh" -n $n -t "$T" -o "subfam_n$n" -x "$S" "$S.fasta" > /dev/null
        run "SubFam -n $n" "subfam_n$n/$S.cons.fasta" "subfam_n$n/$S.chunks.tsv" "$t"
    done

    for id in 0.95 0.90 0.80; do
        # VSEARCH: greedy centroid clustering, with its built-in cluster consensus
        t=$(now)
        vsearch --cluster_fast "$S.fasta" --id $id --threads "$T" --quiet \
            --uc "$S.vs$id.uc" --consout "$S.vs$id.cons.fa"
        awk -F'\t' '$1 == "S" { print $9 "\t" $9 } $1 == "H" { print $9 "\t" $10 }' "$S.vs$id.uc" > "$S.vs$id.tsv"
        sed -i 's/^>centroid=\([^;]*\);.*/>\1/' "$S.vs$id.cons.fa"
        run "VSEARCH --id $id --consout" "$S.vs$id.cons.fa" "$S.vs$id.tsv" "$t"

        # CD-HIT-EST: greedy clustering, representatives are real copies
        t=$(now)
        word=$(awk -v c=$id 'BEGIN { print (c >= 0.9) ? 8 : 5 }')
        cd-hit-est -i "$S.fasta" -o "$S.cdhit$id" -c $id -n "$word" -T "$T" -M 0 -d 0 > /dev/null
        awk '/^>Cluster/ { for (m in mem) print m "\t" rep; delete mem; next }
             { id = $3; sub(/^>/, "", id); sub(/\.\.\.$/, "", id); mem[id] = 1; if ($NF == "*") rep = id }
             END { for (m in mem) print m "\t" rep }' "$S.cdhit$id.clstr" > "$S.cdhit$id.tsv"
        run "CD-HIT-EST -c $id (centroids)" "$S.cdhit$id" "$S.cdhit$id.tsv" "$t"
        "$HERE/cluster_consensus.sh" "$S.fasta" "$S.cdhit$id.tsv" "$S.cdhit$id.cons.fa"
        run "CD-HIT-EST -c $id + mafft/cons" "$S.cdhit$id.cons.fa" "$S.cdhit$id.tsv" "$t"

        # MMseqs2: cascaded clustering, representatives are real copies
        t=$(now)
        mmseqs easy-cluster "$S.fasta" "$S.mm$id" "tmp_mm" --min-seq-id $id -c 0.8 --threads "$T" -v 0 > /dev/null
        awk -F'\t' '{ print $2 "\t" $1 }' "${S}.mm${id}_cluster.tsv" > "$S.mm$id.tsv"
        run "MMseqs2 --min-seq-id $id (centroids)" "${S}.mm${id}_rep_seq.fasta" "$S.mm$id.tsv" "$t"
    done

    # Threshold sweep for the strongest competitor: does some --id work for old families?
    for id in 0.85 0.75 0.70 0.65; do
        t=$(now)
        vsearch --cluster_fast "$S.fasta" --id $id --threads "$T" --quiet \
            --uc "$S.vs$id.uc" --consout "$S.vs$id.cons.fa"
        awk -F'\t' '$1 == "S" { print $9 "\t" $9 } $1 == "H" { print $9 "\t" $10 }' "$S.vs$id.uc" > "$S.vs$id.tsv"
        sed -i 's/^>centroid=\([^;]*\);.*/>\1/' "$S.vs$id.cons.fa"
        run "VSEARCH --id $id --consout" "$S.vs$id.cons.fa" "$S.vs$id.tsv" "$t"
    done
    } | sed -u "s/^/$S\t/" | tee -a results.tsv
done
