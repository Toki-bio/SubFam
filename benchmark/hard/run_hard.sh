#!/usr/bin/env bash
# Harder simulation (simulate2.py), SubFam variants vs VSEARCH thresholds, several seeds.
# usage: run_hard.sh OUTDIR "SEEDS" [THREADS]      e.g. run_hard.sh hard_out "1 2 3 4" 2
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); SUB="$HERE/../../SubFam.sh"; OUT=$1; SEEDS=$2; T=${3:-2}
mkdir -p "$OUT"; cd "$OUT"
# scenario: name | simulate2.py options. mixed = four young and four old subfamilies in one dataset.
SCEN=("young|--divergence 0.03" "middle|--divergence 0.08" "old|--divergence 0.15"
      "mixed|--ages 0.03,0.15,0.03,0.15,0.03,0.15,0.03,0.15" "heavytrunc|--divergence 0.08 --trunc-frac 0.8")
printf 'seed\tscenario\tmethod\treps\treps_ge10\tsingleton\tpurity\tid_all\tid_noCpG\trec_all\trec_noCpG\n' > results.tsv
for seed in $SEEDS; do for sc in "${SCEN[@]}"; do
  name=${sc%%|*}; opts=${sc#*|}; w=s${seed}_$name; mkdir -p $w; cd $w
  python3 "$HERE/simulate2.py" $name $opts --seed $seed
  ev() { python3 "$HERE/evaluate2.py" "$1" "$2" "$3" $name.truth.tsv $name.masters.fasta | sed "s/^/$seed\t$name\t/" >> ../results.tsv; }
  for v in "-n 50" "-n 20" "-n 20 -c" "-m -n 20"; do
    d=sf${v// /}; "$SUB" $v -t $T -o $d -x $name $name.fasta > /dev/null 2>&1
    ev "SubFam $v" $d/$name.cons.fasta $d/$name.chunks.tsv
  done
  for id in 0.95 0.90 0.85 0.80 0.75 0.70; do
    vsearch --cluster_fast $name.fasta --id $id --threads $T --quiet --uc v.uc --consout v.cons.fa
    awk -F'\t' '$1=="S"{print $9"\t"$9} $1=="H"{print $9"\t"$10}' v.uc > v.tsv; sed -i 's/^>centroid=\([^;]*\);.*/>\1/' v.cons.fa
    ev "VSEARCH $id" v.cons.fa v.tsv
  done
  cd ..
done; done
echo finished > done.flag
