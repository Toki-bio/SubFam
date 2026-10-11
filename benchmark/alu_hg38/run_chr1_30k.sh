#!/usr/bin/env bash
# hg38 chr1 Alu stress case: 30,000 random full-length copies (seed 1); SubFam -n 50 and -n 20 (k-mer order), -n 50 -P; wall time and peak RSS
set -e; SUB=/home/user/SubFam/SubFam.sh; cd "$(dirname "$0")"
printf 'run\tsec\tpeak_rss_MB\trows\n' > runs.tsv
for v in "n50:-n 50" "n20:-n 20" "n50P:-n 50 -P"; do
  tag=${v%%:*}; opt=${v#*:}
  python3 - "$SUB" $opt -t 8 -o $tag -x alu s1.fa > $tag.time <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=open(sys.argv[1+sys.argv.index('-o')]+'.stdout','w') if False else subprocess.DEVNULL, stderr=subprocess.STDOUT)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d %d' % (time.time() - t, ru.ru_maxrss // 1024, r.returncode))
PY
  read sec mb rc < $tag.time; rows=$(grep -c '>' $tag/alu.cons.fasta 2>/dev/null || echo 0)
  printf '%s\t%s\t%s\t%s\n' "$tag" "$sec" "$mb" "$rows" >> runs.tsv
done
echo done > done.flag
