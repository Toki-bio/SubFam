#!/usr/bin/env bash
# Noise floor for the flank test: element-only copies in shuffled input order; and Vassetzky 2021 consensuses as references.
set -euo pipefail
W=/staging/tmp/subfam_flanks
cd "$W" || exit 1
export TMPDIR=$W/tmp
TOP=/staging/tmp/subfam_b2_purity/SubFam; T=2
cp "$TOP/benchmark/sine_b2/kramerov2021_B2_consensus.fa" .
for s in 1 2; do
    python3 - f0.fa s$s.fa $s <<'PY'
import sys, random
src, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs = open(src).read().split('>')[1:]; random.Random(seed).shuffle(recs)
open(out, 'w').write(''.join('>' + r for r in recs))
PY
    /usr/bin/time -f "%e %M" -o s$s.time "$TOP/SubFam.sh" -n 50 -t $T -o sf_s$s -x s$s s$s.fa > sf_s$s.log 2>&1
done
python3 - <<'PY'
import importlib.util, sys
src = open('compare.py').read().split("p0 = part(")[0]          # reuse the functions only
exec(src)
p0 = part('sf_f0/f0.chunks.tsv'); r0 = fa('sf_f0/f0.cons.fasta')
refs = {'Dfam': fa('dfam_B2.fa'), 'V2021': fa('kramerov2021_B2_consensus.fa')}
print('run\tpairs_kept_vs_f0\tARI_vs_f0\tf0_rows_found_>=0.95\tmedian\tmin\tDfam_found\tDfam_ids\tV2021_found\tV2021_ids')
for run in ('f0', 's1', 's2', 'f25', 'f100', 'f500'):
    p = part('sf_%s/%s.chunks.tsv' % (run, run)); r = fa('sf_%s/%s.cons.fasta' % (run, run))
    b0 = sorted(best(s, r) for s in r0.values())
    out = [run, '%.3f' % kept_pairs(p0, p), '%.3f' % ari(p0, p), '%d/60' % sum(x >= 0.95 for x in b0), '%.3f' % b0[30], '%.3f' % b0[0]]
    for k in ('Dfam', 'V2021'):
        bd = [best(s, r) for s in refs[k].values()]
        out += ['%d/%d' % (sum(x >= 0.95 for x in bd), len(bd)), ' '.join('%s:%.3f' % (n.split('_')[-1], x) for n, x in zip(refs[k], bd))]
    print('\t'.join(out))
PY
echo ALL_DONE
