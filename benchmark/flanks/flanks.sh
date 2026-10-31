#!/usr/bin/env bash
# Does SubFam's result stand when the copies carry their own unique genomic flanks?
# The same 3,000 mm39 chr19 B2-family copies (make_b2.sh sample, seed 1), cut with 0, 25, 100, 500 bp of flank on both sides.
# Label-free comparison with the flank-0 run: chunk partition agreement and whether the flank-0 rows reappear in the flanked rows.
# DRAGEN rule: at most 4 threads.
set -euo pipefail
W=/staging/tmp/subfam_flanks
mkdir -p "$W"
cd "$W" || exit 1
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
SRC=/staging/tmp/subfam_b2_purity; TOP=$SRC/SubFam; T=4
git -C "$TOP" log -1 --format='SubFam commit %h'
cat > cut.py <<'EOF'
import gzip, sys
rm, fa, ref, F, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
want = {}; order = []; name = None
for l in open(ref):
    if l.startswith('>'): name = l[1:].split()[0]; order.append(name); want[name] = []
    else: want[name].append(l.strip())
want = {k: ''.join(v) for k, v in want.items()}
seq = ''.join(l.strip() for l in gzip.open(fa, 'rt') if not l.startswith('>')).upper()
strand = {}
for l in gzip.open(rm, 'rt'):
    c = l.rstrip('\n').split('\t')
    if c[5] != 'chr19': continue
    k = '%s_chr19_%s_%s' % (c[10], c[6], c[7])
    if k in want: strand[k] = c[9]
comp = str.maketrans('ACGTN', 'TGCAN'); bad = 0
with open(out, 'w') as o:
    for k in order:
        s, e = int(k.split('_')[-2]), int(k.split('_')[-1])
        a, b = max(0, s - F), min(len(seq), e + F)
        t = seq[a:b]
        if strand[k] == '-': t = t.translate(comp)[::-1]; lf, rf = b - e, s - a
        else: lf, rf = s - a, b - e
        inner = t[lf:len(t) - rf]
        if inner != want[k]: bad += 1
        o.write('>%s\n%s\n' % (k, t))
print('flank', F, 'copies', len(order), 'inner part differs from the flank-0 copy:', bad, file=sys.stderr)
if bad: sys.exit(1)
EOF
for F in 0 25 100 500; do
    python3 cut.py "$SRC/dl/rmsk.txt.gz" "$SRC/dl/chr19.fa.gz" "$SRC/b2_chr19_3000.fa" $F f$F.fa
    /usr/bin/time -f "%e %M" -o f$F.time "$TOP/SubFam.sh" -n 50 -t $T -o sf_f$F -x f$F f$F.fa > sf_f$F.log 2>&1
done
cp "$SRC/dfam_B2.fa" .
cat > compare.py <<'EOF'
import collections, itertools, math
from Bio import Align
al = Align.PairwiseAligner(); al.mode = 'local'; al.match_score = 2; al.mismatch_score = -3; al.open_gap_score = -5; al.extend_gap_score = -2
def fa(f):
    d = collections.OrderedDict(); n = None
    for l in open(f):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v).upper().replace('-', '') for k, v in d.items()}
def part(f): return {l.split()[0]: l.split()[1] for l in open(f)}
def ari(p, q):
    ks = list(p); n = len(ks); c = collections.Counter((p[k], q[k]) for k in ks)
    a = collections.Counter(p[k] for k in ks); b = collections.Counter(q[k] for k in ks); C = lambda x: x * (x - 1) / 2
    sij = sum(C(v) for v in c.values()); sa = sum(C(v) for v in a.values()); sb = sum(C(v) for v in b.values()); e = sa * sb / C(n)
    return (sij - e) / ((sa + sb) / 2 - e)
def kept_pairs(p, q):
    g = collections.defaultdict(list)
    for k, v in p.items(): g[v].append(k)
    tot = same = 0
    for m in g.values():
        for x, y in itertools.combinations(m, 2): tot += 1; same += q[x] == q[y]
    return same / tot
def ident(qs, ts):
    """identity of the best local alignment of query qs inside target ts, over the aligned columns; coverage of qs"""
    a = al.align(qs, ts)[0]; cols = m = 0; qcov = 0
    for (q0, q1), (t0, t1) in zip(*a.aligned):
        m += sum(x == y and x != 'N' for x, y in zip(qs[q0:q1], ts[t0:t1])); qcov += q1 - q0
    s = a.aligned; span_q = s[0][-1][1] - s[0][0][0]; span_t = s[1][-1][1] - s[1][0][0]
    cols = max(span_q, span_t)
    return m / cols, qcov / len(qs)
def best(qs, rows):
    b = (0, 0)
    for r in rows.values():
        i, c = ident(qs, r)
        if c >= 0.9 and i > b[0]: b = (i, c)
    return b[0]
p0 = part('sf_f0/f0.chunks.tsv'); r0 = fa('sf_f0/f0.cons.fasta'); dfam = fa('dfam_B2.fa')
print('flank\twall_s\tpeak_MB\trows\tARI_vs_0\tpairs_kept_vs_0\tf0_rows_found_>=0.95\tf0_rows_median_best_id\tf0_rows_min_best_id\tDfam_found_>=0.95\tDfam_best_ids')
for F in (0, 25, 100, 500):
    p = part('sf_f%d/f%d.chunks.tsv' % (F, F)); r = fa('sf_f%d/f%d.cons.fasta' % (F, F))
    t = open('f%d.time' % F).read().split()
    b0 = sorted(best(s, r) for s in r0.values()); bd = [best(s, r) for s in dfam.values()]
    print('%d\t%s\t%d\t%d\t%.3f\t%.3f\t%d/%d\t%.3f\t%.3f\t%d/5\t%s' % (F, t[0], int(t[1]) // 1024, len(r), ari(p0, p), kept_pairs(p0, p),
          sum(x >= 0.95 for x in b0), len(b0), b0[len(b0) // 2], b0[0], sum(x >= 0.95 for x in bd), ' '.join('%.3f' % x for x in bd)))
EOF
python3 compare.py | tee flanks.tsv
echo ALL_DONE
