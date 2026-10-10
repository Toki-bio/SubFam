#!/usr/bin/env python3
"""vsearch_sweep.py COPIES.fa CHUNKS.tsv TRUTH1.tsv NAME1 TRUTH2.tsv NAME2 [--threads T --ids 0.80,0.90,0.95,0.98]

The threshold baseline: vsearch --cluster_fast (global identity, --iddef 2, plus strand) at each identity, against
SubFam's chunks (CHUNKS.tsv) and one family-wide group. Per method: number of groups (= rows an aligner would get),
groups of one copy, and the weighted purity of the groups for each labelling (fraction of copies that carry the
majority label of their own group; 1.0 = every group is one class). Purity rewards many small groups, so the number
of groups is printed next to it; SubFam's N=50 chunks are not tuned to the labels."""
import collections, subprocess, sys, tempfile, os
args = sys.argv[1:]; threads, ids = 4, '0.60,0.70,0.80,0.90,0.95,0.98'
for f in ('--threads', '--ids'):
    if f in args: i = args.index(f); v = args[i + 1]; del args[i:i + 2]; threads, ids = (int(v), ids) if f == '--threads' else (threads, v)
copies, chunks, t1, n1, t2, n2 = args[:6]
truth = [dict(l.split()[:2] for l in open(t)) for t in (t1, t2)]
names = [l[1:].split()[0] for l in open(copies) if l.startswith('>')]
BIG = 10                                             # a group of at least this many copies is a 'real' group
def score(groups):                                   # groups: id -> list of copies
    out = []
    for tr in truth:
        pur = sum(collections.Counter(tr[m] for m in g).most_common(1)[0][1] for g in groups.values())
        out.append(pur / sum(len(g) for g in groups.values()))
    return out
def row(label, groups):
    p = score(groups); sizes = sorted(len(g) for g in groups.values())
    big = {k: g for k, g in groups.items() if len(g) >= BIG}
    cov = sum(len(g) for g in big.values()) / len(names); pb = score(big) if big else [float('nan')] * 2
    print('%-26s %6d %6d %7d | %5.1f%% %4d %.3f %.3f | %.3f %.3f' % (label, len(groups), sizes[len(sizes) // 2], sum(s == 1 for s in sizes), 100 * cov, len(big), pb[0], pb[1], p[0], p[1]))
print('%-26s %6s %6s %7s | %s (groups of >= %d copies) | all groups' % ('method', 'groups', 'med.sz', 'single', 'copies covered, n groups, purity %s, %s' % (n1, n2), BIG))
row('one family-wide consensus', {'all': names})
g = collections.defaultdict(list)
for l in open(chunks):
    f = l.split()
    if len(f) >= 2: g[f[1]].append(f[0])
row('SubFam chunks (%d)' % len(g), g)
with tempfile.TemporaryDirectory() as t:
    for i in ids.split(','):
        subprocess.run(['vsearch', '--cluster_fast', copies, '--id', i, '--iddef', '2', '--strand', 'plus', '--uc', t + '/c.uc',
                        '--threads', str(threads), '--quiet'], check=True)
        g = collections.defaultdict(list)
        for l in open(t + '/c.uc'):
            f = l.split('\t')
            if f[0] in 'SH': g[f[1]].append(f[8].split()[0])
        row('vsearch --id %s' % i, g)
