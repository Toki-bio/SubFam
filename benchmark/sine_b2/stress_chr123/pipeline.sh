#!/usr/bin/env bash
# B2 stress case, mm39 chr1-3: download, extract full-length B2-family copies, sample 30,000 (seed 1),
# then SubFam -n 50 (k-mer), -n 50 -P, -n 20; wall time and peak RSS of each run.
#   pipeline.sh [WORKDIR] [THREADS]     (WORKDIR default: this directory, where the committed outputs were made)
# The 2026-10-11 run shared the CPU with the Alu runs for -n 50 (13 min); on an idle machine it takes about 6.3 min (seeds.sh).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); R=$(cd "$HERE/../.." && pwd); SUB=$R/../SubFam.sh
W=${1:-$HERE}; T=${2:-8}; mkdir -p "$W/dl"; cd "$W"
U=https://hgdownload.soe.ucsc.edu/goldenPath/mm39
[ -s dl/rmsk.txt.gz ] || curl -sSf -o dl/rmsk.txt.gz $U/database/rmsk.txt.gz
for c in chr1 chr2 chr3; do
  [ -s dl/$c.fa.gz ] || curl -sSf -o dl/$c.fa.gz $U/chromosomes/$c.fa.gz
  python3 -I "$R/line_real/extract_l1.py" dl/rmsk.txt.gz dl/$c.fa.gz b2_$c.fa 170 '^(B2_Mm1a|B2_Mm1t|B2_Mm2|B3|B3A)$' 2>> extract.log
done
cat b2_chr1.fa b2_chr2.fa b2_chr3.fa > b2_chr123.fa; grep -c '>' b2_chr123.fa > n_all.txt     # 41,989 in the 2026-10-11 run
python3 "$R/sine_b2/subsample.py" b2_chr123.fa 30000 1 s1.fa s1.truth.tsv
printf 'run\tsec\tpeak_rss_MB\trows\n' > runs.tsv
for v in "n50:-n 50" "n50P:-n 50 -P" "n20:-n 20"; do
  tag=${v%%:*}; opt=${v#*:}
  # shellcheck disable=SC2086
  python3 - "$SUB" $opt -t "$T" -o $tag -x b2 s1.fa > $tag.time <<'PY'
import resource, subprocess, sys, time
t = time.time(); r = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
ru = resource.getrusage(resource.RUSAGE_CHILDREN); print('%.0f %d %d' % (time.time() - t, ru.ru_maxrss // 1024, r.returncode))
PY
  read -r sec mb rc < $tag.time; [ "$rc" = 0 ] || { echo "pipeline.sh: SubFam $opt failed (exit $rc)" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\n' "$tag" "$sec" "$mb" "$(grep -c '>' $tag/b2.cons.fasta)" >> runs.tsv
done
