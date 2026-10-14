#!/usr/bin/env bash
# B2 stress case, mm39 chr1-3: wait for the chromosome downloads, extract full-length B2 copies, sample 30,000 (seed 1),
# wait until the Alu runs have released the CPU, then SubFam -n 50 (k-mer), -n 50 -P, -n 20; time and peak RSS.
set -e; cd "$(dirname "$0")"; R=$(cd "$(dirname "$0")/../../.." && pwd)/benchmark; SUB=$R/../SubFam.sh
until [ -f dl/dl.flag ]; do sleep 20; done
for c in chr1 chr2 chr3; do python3 -I $R/line_real/extract_l1.py ../b2/dl/rmsk.txt.gz dl/$c.fa.gz b2_$c.fa 170 '^(B2_Mm1a|B2_Mm1t|B2_Mm2|B3|B3A)$' 2>> extract.log; done
cat b2_chr1.fa b2_chr2.fa b2_chr3.fa > b2_chr123.fa; grep -c '>' b2_chr123.fa > n_all.txt
python3 $R/sine_b2/subsample.py b2_chr123.fa 30000 1 s1.fa s1.truth.tsv
until [ -f ../alu_hg38/done.flag ]; do sleep 30; done
printf 'run\tsec\tpeak_rss_MB\trows\n' > runs.tsv
for v in "n50:-n 50" "n50P:-n 50 -P" "n20:-n 20"; do
  tag=${v%%:*}; opt=${v#*:}
  python3 - "$SUB" $opt -t 8 -o $tag -x b2 s1.fa > $tag.time <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d %d' % (time.time() - t, ru.ru_maxrss // 1024, r.returncode))
PY
  read sec mb rc < $tag.time; rows=$(grep -c '>' $tag/b2.cons.fasta 2>/dev/null || echo 0)
  printf '%s\t%s\t%s\t%s\n' "$tag" "$sec" "$mb" "$rows" >> runs.tsv
done
echo done > done.flag
