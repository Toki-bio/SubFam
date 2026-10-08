#!/usr/bin/env python3
"""Score representatives of a LINE-like simulation (long masters, partial coverage allowed).

usage: evaluate_line.py METHOD REPS.fasta MEMBERS.tsv TRUTH.tsv MASTERS.fasta
Each representative is aligned to the masters with vsearch (--usearch_global, global alignment, terminal gaps not counted); identity is taken over the aligned span,
leaving out the master's CpG positions (CpG decay moves a consensus off the master there). Coverage = master positions between the first and last column that carries a query base / master length (vsearch's tilo/tihi include terminal gaps of a global alignment and would give 1.0).
A subfamily is recovered at coverage C when a representative built from >= 3 copies, whose group is mostly that subfamily, reaches id_noCpG >= 99 %
and covers >= C of the master (C = 0.2, 0.5, 0.9: about the 3' 1.2 kb, the 3' half, nearly all).
prints: method reps reps_ge10 singleton purity median_cov_ge10 median_idnoCpG_ge10 rec@0.2 rec@0.5 rec@0.9 (each of 8)"""
import collections, statistics, subprocess, sys, tempfile

def read_fa(p):
    d = {}; cur = None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'): cur = l[1:].split()[0]; d[cur] = []
        else: d[cur].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}

method, reps_f, members, truth_f, masters_f = sys.argv[1:6]
truth = dict(l.split() for l in open(truth_f))
groups = collections.defaultdict(list)
for line in open(members):
    sid, gid = line.split('\t')[:2]
    groups[gid.strip()].append(truth[sid])
assert sum(map(len, groups.values())) == len(truth)
majority = {g: collections.Counter(v).most_common(1)[0] for g, v in groups.items()}
purity = sum(n for _, n in majority.values()) / len(truth)
single = sum(len(v) for v in groups.values() if len(v) == 1) / len(truth)
masters = read_fa(masters_f); reps = read_fa(reps_f)
big = [r for r in reps if len(groups[r]) >= 10]
cpg = {k: {j for i in range(len(m) - 1) if m[i:i + 2] == 'CG' for j in (i, i + 1)} for k, m in masters.items()}
with tempfile.NamedTemporaryFile('r', suffix='.tsv') as h:
    subprocess.run(['vsearch', '--usearch_global', reps_f, '--db', masters_f, '--id', '0.5', '--iddef', '2', '--strand', 'both',
                    '--maxaccepts', '0', '--maxrejects', '0', '--maxhits', '1', '--userout', h.name,
                    '--userfields', 'query+target+id+tilo+tihi+qrow+trow', '--quiet', '--threads', '2'], check=True)
    res = {}
    for l in h:
        q, t, ident, lo, hi, qrow, trow = l.rstrip('\n').split('\t'); qrow, trow = qrow.upper(), trow.upper()
        pos = int(lo) - 1; same = diff = 0
        first = next((i for i, c in enumerate(qrow) if c != '-'), 0)
        last = max((i for i, c in enumerate(qrow) if c != '-'), default=len(qrow) - 1)
        for i, (qc, tc) in enumerate(zip(qrow, trow)):
            if tc != '-':
                if first <= i <= last and pos not in cpg[t]:
                    if qc == tc: same += 1
                    else: diff += 1
                pos += 1
            elif first <= i <= last:
                diff += 1
        cov_cols = sum(1 for i, (qc, tc) in enumerate(zip(qrow, trow)) if first <= i <= last and tc != '-')   # master positions between the first and last aligned column
        res[q] = (t, 100.0 * same / max(1, same + diff), cov_cols / len(masters[t]))
def rec(c):
    return len({res[r][0] for r in res if len(groups[r]) >= 3 and res[r][1] >= 99 and res[r][2] >= c and majority[r][0] == res[r][0]})
bigres = [res[r] for r in big if r in res]
print('\t'.join([method, str(len(reps)), str(len(big)), '%.2f' % single, '%.3f' % purity,
                 '%.2f' % (statistics.median(x[2] for x in bigres) if bigres else 0), '%.1f' % (statistics.median(x[1] for x in bigres) if bigres else 0),
                 '%d/8' % rec(0.2), '%d/8' % rec(0.5), '%d/8' % rec(0.9)]))
