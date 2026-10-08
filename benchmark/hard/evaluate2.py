#!/usr/bin/env python3
"""Score representatives against the simulated truth, CpG-aware and coverage-aware.

usage: evaluate2.py METHOD REPS.fasta MEMBERS.tsv TRUTH.tsv MASTERS.fasta [COV=0.95]

Identity of a representative to a true subfamily master is computed from a local alignment (match 2, mismatch -3, gap -5/-2, both strands),
  id_all     over all aligned columns, as evaluate.py does
  id_noCpG   leaving out the master's CpG positions: CpG decay is expected to move every consensus away from the master there
             (the consensus of decayed copies is TG/CA), so only non-CpG columns say whether the row is the master.
Coverage = aligned span of the master / master length.
A subfamily is recovered when a representative whose group is mostly that subfamily reaches id >= 99 % and coverage >= COV.
prints: method reps reps_ge10 singleton_frac purity median_id_all_ge10 median_id_noCpG_ge10 recovered(id_all) recovered(id_noCpG)  (each of 8, with coverage)"""
import collections, statistics, sys
from Bio import Align
from Bio.Seq import Seq

def read_fa(p):
    d = {}; cur = None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'): cur = l[1:].split()[0]; d[cur] = []
        else: d[cur].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}

method, reps_f, members, truth_f, masters_f = sys.argv[1:6]
cov_min = float(sys.argv[6]) if len(sys.argv) > 6 else 0.95
truth = dict(l.split() for l in open(truth_f))
groups = collections.defaultdict(list)
for line in open(members):
    sid, gid = line.split('\t')[:2]
    groups[gid.strip()].append(truth[sid])
assert sum(map(len, groups.values())) == len(truth)
majority = {g: collections.Counter(v).most_common(1)[0] for g, v in groups.items()}
purity = sum(n for _, n in majority.values()) / len(truth)
single = sum(len(v) for v in groups.values() if len(v) == 1) / len(truth)
reps = read_fa(reps_f); masters = read_fa(masters_f)
big = [r for r in reps if len(groups[r]) >= 10]
cpg = {}
for k, m in masters.items():
    s = set()
    for i in range(len(m) - 1):
        if m[i:i + 2] == 'CG': s.update((i, i + 1))
    cpg[k] = s
al = Align.PairwiseAligner(); al.mode = 'local'
al.match_score = 2; al.mismatch_score = -3; al.open_gap_score = -5; al.extend_gap_score = -2

def stats(rep, name):
    m = masters[name]; best = None; bq = None
    for q in (rep, str(Seq(rep).reverse_complement())):
        a = al.align(m, q)[0]
        if best is None or a.score > best.score: best, bq = a, q
    tb, qb = best.aligned
    same = diff = same_nc = diff_nc = gap = 0
    for (t0, t1), (q0, q1) in zip(tb, qb):
        for k in range(t1 - t0):
            eq = m[t0 + k] == bq[q0 + k]
            same += eq; diff += (not eq)
            if (t0 + k) not in cpg[name]:
                same_nc += eq; diff_nc += (not eq)
    for (a0, a1), (b0, b1) in zip(zip(tb[:-1], qb[:-1]), zip(tb[1:], qb[1:])):
        gap += max(b0[0] - a0[1], b1[0] - a1[1])
    id_all = 100.0 * same / max(1, same + diff + gap)
    id_nc = 100.0 * same_nc / max(1, same_nc + diff_nc + gap)
    cov = (tb[-1][1] - tb[0][0]) / len(m)
    return id_all, id_nc, cov

best = {}
for r, s in reps.items():
    if len(groups[r]) < 3:      # a representative built from fewer than 3 copies cannot be a recovered consensus (saves time)
        continue
    cand = [(name,) + stats(s, name) for name in masters]
    best[r] = {'all': max(cand, key=lambda x: x[1]), 'nc': max(cand, key=lambda x: x[2])}
def med(key, idx):
    v = [best[r][key][idx] for r in big if r in best]
    return statistics.median(v) if v else float('nan')
rec_all = {best[r]['all'][0] for r in best if best[r]['all'][1] >= 99 and best[r]['all'][3] >= cov_min and majority[r][0] == best[r]['all'][0]}
rec_nc = {best[r]['nc'][0] for r in best if best[r]['nc'][2] >= 99 and best[r]['nc'][3] >= cov_min and majority[r][0] == best[r]['nc'][0]}
print('\t'.join([method, str(len(reps)), str(len(big)), '%.2f' % single, '%.3f' % purity, '%.1f' % med('all', 1), '%.1f' % med('nc', 2), '%d/8' % len(rec_all), '%d/8' % len(rec_nc)]))
