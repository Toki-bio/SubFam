#!/usr/bin/env bash
B=/home/user/SubFam/benchmark; S=/home/user/SubFam/SubFam.sh
printf 'seed\tscenario\tmethod\treps\treps_ge10\tsingleton\tpurity\tmedian_id\trecovered\n' > results.tsv
for seed in 1 2 3 4 5 6 7 8; do
 for scen in young:0.03 middle:0.08 old:0.15; do
  s=${scen%%:*}; d=${scen#*:}; w=s${seed}_$s; mkdir -p $w; cd $w
  python3 $B/simulate.py $s --divergence $d --seed $seed
  for n in 50 20; do
    $S -n $n -t 4 -o sf$n -x $s $s.fasta >/dev/null 2>&1
    python3 $B/evaluate.py "SubFam -n $n" sf$n/$s.cons.fasta sf$n/$s.chunks.tsv $s.truth.tsv $s.masters.fasta | sed "s/^/$seed\t$s\t/" >> ../results.tsv
  done
  for id in 0.95 0.90 0.85 0.80 0.75 0.70 0.65; do
    vsearch --cluster_fast $s.fasta --id $id --threads 4 --quiet --uc v.uc --consout v.cons.fa
    awk -F'\t' '$1=="S"{print $9"\t"$9} $1=="H"{print $9"\t"$10}' v.uc > v.tsv; sed -i 's/^>centroid=\([^;]*\);.*/>\1/' v.cons.fa
    python3 $B/evaluate.py "VSEARCH $id" v.cons.fa v.tsv $s.truth.tsv $s.masters.fasta | sed "s/^/$seed\t$s\t/" >> ../results.tsv
  done
  cd ..
 done
done
echo finished > done.flag
