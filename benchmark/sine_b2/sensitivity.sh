#!/usr/bin/env bash
# usage: sensitivity.sh B2_WORKDIR   (the directory of make_b2.sh: b2_chr19.fa, b2_chr19_3000.fa, b2_truth.tsv, b2_truth_kramerov.tsv); output sens.tsv, committed as sensitivity.tsv
# sensitivity of the B2 result to N, k, p and to the sample: purity of chunks vs rmsk and vs the 2021 labels
set -e; B=$1; SUB=/home/user/SubFam/SubFam.sh
printf 'setting\tsample\tchunks\tpurity_rmsk\tpurity_v2021\tmean_agree\n' > sens.tsv
score() { # dir prefix sample setting
  python3 - "$1/$2.chunks.tsv" "$B/b2_truth.tsv" "$3" "$1/$2.chunkstats.tsv" "$4" "$5" <<'PY'
import sys, collections
ch, tr, tk, st, setting, sample = sys.argv[1:7]
m = collections.defaultdict(list)
for l in open(ch): f = l.split(); m[f[1]].append(f[0])
def pur(tf):
    t = dict(l.split()[:2] for l in open(tf)); return sum(collections.Counter(t[x] for x in g).most_common(1)[0][1] for g in m.values()) / sum(len(g) for g in m.values())
ag = [float(l.split('\t')[5]) for l in open(st) if not l.startswith('chunk')]
print('%s\t%s\t%d\t%.3f\t%.3f\t%.3f' % (setting, sample, len(m), pur(tr), pur(tk), sum(ag) / len(ag)))
PY
}
for n in 20 30 50 100; do $SUB -n $n -t 8 -o r_n$n -x b2 $B/b2_chr19_3000.fa >/dev/null; score r_n$n b2 $B/b2_truth_kramerov.tsv "n=$n k=6 p=0.36" seed1 >> sens.tsv; done
for k in 4 5 8 10; do $SUB -n 50 -k $k -t 8 -o r_k$k -x b2 $B/b2_chr19_3000.fa >/dev/null; score r_k$k b2 $B/b2_truth_kramerov.tsv "n=50 k=$k p=0.36" seed1 >> sens.tsv; done
for p in 0.25 0.5; do $SUB -n 50 -p $p -t 8 -o r_p$p -x b2 $B/b2_chr19_3000.fa >/dev/null; score r_p$p b2 $B/b2_truth_kramerov.tsv "n=50 k=6 p=$p" seed1 >> sens.tsv; done
$SUB -n 50 -m -t 8 -o r_m -x b2 $B/b2_chr19_3000.fa >/dev/null; score r_m b2 $B/b2_truth_kramerov.tsv "n=50 mafft-order p=0.36" seed1 >> sens.tsv
# other samples of the 5,972 copies
for seed in 2 3 4; do
  python3 /home/user/SubFam/benchmark/sine_b2/subsample.py $B/b2_chr19.fa 3000 $seed s$seed.fa s$seed.truth.tsv
  python3 /home/user/SubFam/benchmark/sine_b2/label_by_consensus.py s$seed.fa /home/user/SubFam/benchmark/sine_b2/kramerov2021_B2_consensus.fa s$seed.kram.tsv --threads 8 2>/dev/null
  $SUB -n 50 -t 8 -o r_s$seed -x b2 s$seed.fa >/dev/null
  python3 - r_s$seed/b2.chunks.tsv s$seed.truth.tsv s$seed.kram.tsv r_s$seed/b2.chunkstats.tsv "n=50 k=6 p=0.36" seed$seed <<'PY' >> sens.tsv
import sys, collections
ch, tr, tk, st, setting, sample = sys.argv[1:7]
m = collections.defaultdict(list)
for l in open(ch): f = l.split(); m[f[1]].append(f[0])
def pur(tf):
    t = dict(l.split()[:2] for l in open(tf)); return sum(collections.Counter(t[x] for x in g).most_common(1)[0][1] for g in m.values()) / sum(len(g) for g in m.values())
ag = [float(l.split('\t')[5]) for l in open(st) if not l.startswith('chunk')]
print('%s\t%s\t%d\t%.3f\t%.3f\t%.3f' % (setting, sample, len(m), pur(tr), pur(tk), sum(ag) / len(ag)))
PY
done
cat sens.tsv
