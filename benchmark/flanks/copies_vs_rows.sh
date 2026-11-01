#!/usr/bin/env bash
# Copies vs the rows built from them: identity of each to the nearest published consensus (Dfam, Vassetzky 2021). Needs the flank-0 run of flanks.sh. 1 thread.
set -euo pipefail
cd /staging/tmp/subfam_flanks || exit 1
python3 - <<'PY'
import statistics as st
src = open('compare.py').read().split("p0 = part(")[0]
exec(src)
copies = fa('f0.fa'); rows = fa('sf_f0/f0.cons.fasta')
for name, f in (('Dfam', 'dfam_B2.fa'), ('V2021', 'kramerov2021_B2_consensus.fa')):
    refs = fa(f)
    def near(s):
        b = 0
        for r in refs.values():
            i, c = ident(s, r)
            if c >= 0.5 and i > b: b = i
        return b
    ci = [near(s) for s in copies.values()]; ri = [near(s) for s in rows.values()]; print(name, 'copies with no hit', sum(x == 0 for x in ci), 'rows with no hit', sum(x == 0 for x in ri)); ci = [x for x in ci if x]; ri = [x for x in ri if x]
    q = lambda v: '%.3f %.3f %.3f' % (st.quantiles(v, n=4)[0], st.median(v), st.quantiles(v, n=4)[2])
    print(name, 'copies n=%d q1/med/q3 %s  mean %.3f | rows n=%d q1/med/q3 %s mean %.3f' % (len(ci), q(ci), st.mean(ci), len(ri), q(ri), st.mean(ri)))
PY
