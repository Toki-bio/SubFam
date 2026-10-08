#!/usr/bin/env python3
"""Consensus recovery: does a set of chunk consensuses contain each reference consensus?

usage: recover.py REFS.fa CONS.fa OUT_PREFIX [MIN_COV]
  REFS.fa  reference consensuses (Dfam Alu families, Price 213 ...)
  CONS.fa  chunk consensuses (SubFam PREFIX.cons.fasta)
For every reference the best chunk consensus is found by local alignment (match 2, mismatch -3,
gap open -5, extend -2; both orientations tried). Identity = identical columns / alignment columns
(gaps count as columns); coverage = aligned reference span / reference length.
A reference is 'recovered at T' if some chunk consensus aligns with coverage >= MIN_COV (default 0.90)
and identity >= T. T = 0.95, 0.98, 0.99 are reported (fixed in advance, see README).
Also reports, per chunk consensus, its nearest reference (redundancy / novelty).
Writes OUT_PREFIX.ref_best.tsv, OUT_PREFIX.cons_best.tsv and prints a summary."""
import sys, collections
from Bio import Align
from Bio.Seq import Seq

def read_fa(p):
    names, seqs = [], []
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'):
            names.append(l[1:].split()[0]); seqs.append([])
        elif names:
            seqs[-1].append(l.strip())
    return names, [''.join(s).upper().replace('-', '').replace('N', '') for s in seqs]

refs_f, cons_f, out = sys.argv[1:4]
min_cov = float(sys.argv[4]) if len(sys.argv) > 4 else 0.90
rn, rs = read_fa(refs_f); cn, cs = read_fa(cons_f)
al = Align.PairwiseAligner(); al.mode = 'local'
al.match_score = 2; al.mismatch_score = -3; al.open_gap_score = -5; al.extend_gap_score = -2

def best_hit(ref, cons):
    """-> (identity, coverage) of the best local alignment of ref against cons, either strand."""
    best = (0.0, 0.0, -1)
    for q in (cons, str(Seq(cons).reverse_complement())):
        a = al.align(ref, q)[0]
        c = a.counts()
        cols = a.shape[1]
        ident = c.identities / cols if cols else 0.0
        cov = (a.aligned[0][-1][1] - a.aligned[0][0][0]) / len(ref)
        if a.score > best[2]:
            best = (ident, cov, a.score)
    return best[0], best[1]

ref_best = []
cons_best = collections.defaultdict(lambda: ('', 0.0, 0.0))
for i, r in enumerate(rs):
    top = ('', 0.0, 0.0)
    for j, c in enumerate(cs):
        if not c:
            continue
        ident, cov = best_hit(r, c)
        if cov >= min_cov and ident > top[1]:
            top = (cn[j], ident, cov)
        if cov >= min_cov and ident > cons_best[cn[j]][1]:
            cons_best[cn[j]] = (rn[i], ident, cov)
    ref_best.append((rn[i], len(r)) + top)
with open(out + '.ref_best.tsv', 'w') as f:
    f.write('reference\tref_len\tbest_consensus\tidentity\tcoverage\n')
    for r in ref_best:
        f.write('%s\t%d\t%s\t%.4f\t%.3f\n' % r)
with open(out + '.cons_best.tsv', 'w') as f:
    f.write('consensus\tnearest_reference\tidentity\tcoverage\n')
    for c in cn:
        r = cons_best.get(c, ('', 0.0, 0.0))
        f.write('%s\t%s\t%.4f\t%.3f\n' % (c, r[0], r[1], r[2]))
n = len(ref_best)
print('references: %d   chunk consensuses: %d   min coverage %.2f' % (n, len(cn), min_cov))
for t in (0.95, 0.98, 0.99):
    k = sum(1 for r in ref_best if r[3] >= t)
    print('recovered at identity >= %.2f: %d / %d (%.1f %%)' % (t, k, n, 100.0 * k / n if n else 0))
ids = sorted(r[3] for r in ref_best)
print('best identity per reference: median %.3f, min %.3f' % (ids[len(ids) // 2], ids[0]) if ids else '')
