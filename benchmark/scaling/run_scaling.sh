#!/usr/bin/env bash
# Wall time and peak memory of the whole SubFam pipeline against input size, per ordering method.
#   run_scaling.sh FASTA OUTDIR [SIZES="1000 3000 10000 30000 100000"] [THREADS=16]
# FASTA is a real SINE copy set (any size >= the largest SIZE). Each size is a random subsample (seed 11).
# Methods: kmer (C, default), mafft (-m, MAFFT guide tree), parttree (-P), old (SubFam.old.sh, the original
# script; only up to 30,000). Needs: seqkit, mafft, GNU time (/usr/bin/time -v). Output: OUTDIR/scaling.tsv.
set -uo pipefail
FA=${1:?fasta}
OUT=${2:?outdir}
SIZES=${3:-"1000 3000 10000 30000 100000"}
THREADS=${4:-16}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBFAM=${SUBFAM:-$HERE/../../SubFam.sh}
OLD=${OLD:-}
mkdir -p "$OUT" && cd "$OUT" || exit 1
printf 'copies\tmethod\twall_s\tpeak_rss_mb\tchunks\texit\tload_start\n' > scaling.tsv
for n in $SIZES; do
    seqkit sample -n "$n" -s 11 "$FA" 2> /dev/null | seqkit seq -w 0 > "in_$n.fa"
    got=$(grep -c '^>' "in_$n.fa")
    for m in kmer mafft parttree old; do
        case $m in
            kmer)     cmd=(bash "$SUBFAM" -n 50 -t "$THREADS" -o "run_${n}_$m" "in_$n.fa") ;;
            mafft)    cmd=(bash "$SUBFAM" -n 50 -t "$THREADS" -m -o "run_${n}_$m" "in_$n.fa") ;;
            parttree) cmd=(bash "$SUBFAM" -n 50 -t "$THREADS" -P -o "run_${n}_$m" "in_$n.fa") ;;
            old)      [ -n "$OLD" ] && [ "$n" -le 30000 ] || continue
                      mkdir -p "run_${n}_old"; cp "in_$n.fa" "run_${n}_old/in.fa"
                      cmd=(bash -c "cd run_${n}_old && bash $OLD in.fa 50") ;;
        esac
        load=$(cut -d' ' -f1 /proc/loadavg)
        { /usr/bin/time -f '%e\t%M\t%x' "${cmd[@]}" > "log_${n}_$m.txt" 2> "time_${n}_$m.txt"; } 2>> "time_${n}_$m.txt"
        t=$(tail -1 "time_${n}_$m.txt")
        wall=$(echo "$t" | cut -f1); rss=$(echo "$t" | awk -F'\t' '{printf "%d", $2/1024}'); rc=$(echo "$t" | cut -f3)
        ch=$(ls "run_${n}_$m"/*.chunks.tsv 2> /dev/null | head -1)
        nch=0; [ -n "$ch" ] && nch=$(cut -f2 "$ch" | sort -u | wc -l)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$got" "$m" "$wall" "$rss" "$nch" "$rc" "$load" >> scaling.tsv
        rm -rf "run_${n}_$m"
    done
    rm -f "in_$n.fa"
done
cat scaling.tsv
