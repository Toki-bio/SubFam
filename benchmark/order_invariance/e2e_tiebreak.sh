#!/usr/bin/env bash
# End to end: SubFam with the content-based tie-break, clean input and 3 shuffles. 2 threads.
set -euo pipefail
cd /staging/tmp/subfam_order || exit 1
export TMPDIR=$PWD
rm -rf SubFam_tb; cp -r /staging/tmp/subfam_b2_purity/SubFam SubFam_tb; cp kmer_order_tb.c SubFam_tb/kmer_order.c
for f in b2_chr19_3000 p1 p2 p3; do ./SubFam_tb/SubFam.sh -n 50 -t 2 -o tbo_$f -x $f $f.fa > tbo_$f.log 2>&1; done
python3 - <<'PY'
src = open('/staging/tmp/subfam_flanks/compare.py').read().split("p0 = part(")[0]
exec(src)
P = lambda r: part('tbo_%s/%s.chunks.tsv' % (r, r)); C = lambda r: fa('tbo_%s/%s.cons.fasta' % (r, r))
r0 = C('b2_chr19_3000')
for r in ('p1', 'p2', 'p3'):
    b = sorted(best(s, C(r)) for s in r0.values())
    print(r, 'pairs kept %.3f' % kept_pairs(P('b2_chr19_3000'), P(r)), 'original rows recurring >=0.95: %d/60' % sum(x >= 0.95 for x in b), 'min best id %.3f' % b[0])
PY
echo ALL_DONE
