#!/usr/bin/env bash
# -P (PartTree) vs k-mer order at 10,000 and 20,000 simulated copies: wall time, peak RSS, masters recovered (evaluate.py)
set -e; SUB=/home/user/SubFam/SubFam.sh; B=/home/user/SubFam/benchmark
[ -s mid10k.fasta ] || python3 $B/simulate.py mid10k --divergence 0.08 --seed 7 --sizes 3000,2000,1500,1250,1000,750,300,200
[ -s mid20k.fasta ] || python3 $B/simulate.py mid20k --divergence 0.08 --seed 8 --sizes 6000,4000,3000,2500,2000,1500,600,400
printf 'input\torder\tsec\tpeak_rss_MB\trecovered\trows\n' > scale.tsv
for inp in mid10k mid20k; do for ord in kmer P; do
  o=${inp}_$ord; opt=""; [ $ord = kmer ] || opt="-P"
  python3 - "$SUB" -n 50 $opt -t 8 -o $o -x s $inp.fasta > t.txt <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d' % (time.time() - t, ru.ru_maxrss))
PY
  read sec kb < t.txt
  rec=$(python3 $B/evaluate.py "x" $o/s.cons.fasta $o/s.chunks.tsv $inp.truth.tsv $inp.masters.fasta | tail -1)
  printf '%s\t%s\t%s\t%d\t%s\n' $inp $ord $sec $((kb/1024)) "$rec" >> scale.tsv
done; done
cat scale.tsv
