#!/usr/bin/env bash
# B2 stress case, seeds 2 and 3: new 30,000-copy samples of the 41,989 chr1-3 copies, SubFam -n 50 (k-mer order), timing, scoring.
#   seeds.sh [WORKDIR] [THREADS]     (WORKDIR: the one pipeline.sh used, holding b2_chr123.fa)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); R=$(cd "$HERE/../.." && pwd); SUB=$R/../SubFam.sh
W=${1:-$HERE}; T=${2:-8}; cd "$W"
[ -s b2_chr123.fa ] || { echo "seeds.sh: run pipeline.sh $W first (b2_chr123.fa missing)" >&2; exit 1; }
if [ ! -s dfam_B2.fa ]; then      # Dfam consensuses of the five B2-family entries, as in ../make_b2.sh
  for a in DF000001735 DF000001736 DF000001737 DF000001738 DF000001739; do
    curl -sSf "https://www.dfam.org/api/families/$a/sequence?format=fasta"; echo
  done | grep -v '^$' | sed 's/^>\(DF[0-9.]*\) \(.*\)/>DFAM_\2/' > dfam_B2.fa
fi
for seed in 2 3; do
  python3 "$R/sine_b2/subsample.py" b2_chr123.fa 30000 $seed s$seed.fa s$seed.truth.tsv
  mkdir -p seed$seed; cp s$seed.fa seed$seed/s1.fa; cp s$seed.truth.tsv seed$seed/s1.truth.tsv; cp dfam_B2.fa "$HERE/chr19_rows.fa" seed$seed/
  python3 - "$SUB" -n 50 -t "$T" -o seed$seed/n50 -x b2 seed$seed/s1.fa > seed$seed/n50.time <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d %d' % (time.time() - t, ru.ru_maxrss // 1024, r.returncode))
PY
  read -r _ _ rc < seed$seed/n50.time; [ "$rc" = 0 ] || { echo "seeds.sh: SubFam failed for seed $seed" >&2; exit 1; }
  python3 "$R/sine_b2/score_big.py" "$(pwd)/seed$seed" n50 > seed$seed/score.txt 2>&1
done
