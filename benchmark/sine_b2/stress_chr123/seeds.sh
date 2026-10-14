#!/usr/bin/env bash
# B2 stress case, seeds 2 and 3: new 30,000-copy samples of the 41,989 chr1-3 copies, SubFam -n 50 (k-mer order), timing, scoring
set -e; cd "$(dirname "$0")"; R=$(cd "$(dirname "$0")/../../.." && pwd)/benchmark; SUB=$R/../SubFam.sh
for seed in 2 3; do
  python3 $R/sine_b2/subsample.py b2_chr123.fa 30000 $seed s$seed.fa s$seed.truth.tsv
  mkdir -p seed$seed; cp s$seed.fa seed$seed/s1.fa; cp s$seed.truth.tsv seed$seed/s1.truth.tsv; cp dfam_B2.fa chr19_rows.fa seed$seed/
  python3 - "$SUB" -n 50 -t 8 -o seed$seed/n50 -x b2 seed$seed/s1.fa > seed$seed/n50.time <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d %d' % (time.time() - t, ru.ru_maxrss // 1024, r.returncode))
PY
  python3 $R/sine_b2/score_big.py $(pwd)/seed$seed n50 > seed$seed/score.txt 2>&1
done
echo done > seeds.flag
