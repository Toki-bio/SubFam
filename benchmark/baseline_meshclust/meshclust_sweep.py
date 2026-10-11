#!/usr/bin/env python3
"""meshclust_sweep.py MESHCLUST2_BIN COPIES.fa CHUNKS.tsv TRUTH.tsv NAME --ids 0.6,0.8 [--threads T] [--extra 'args']
MeShClust2 (Girgis lab; CD-HIT-style .clstr output) at several --id values on the same copies, scored exactly as
benchmark/baseline_vsearch/vsearch_sweep.py scores vsearch and SubFam: groups, median size, singletons, copies in groups of >= 10,
purity (fraction of copies carrying the majority label of their group) of those groups and of all groups. Also runs SubFam's chunks
(CHUNKS.tsv) and one family-wide group through the same scorer, so the three rows are comparable. Wall time per run is printed."""
import collections, os, subprocess, sys, tempfile, time
a = sys.argv[1:]; threads, ids, extra = 8, '0.8,0.9', ''
for f in ('--threads', '--ids', '--extra'):
    if f in a: i = a.index(f); v = a[i + 1]; del a[i:i + 2]; threads, ids, extra = (int(v), ids, extra) if f == '--threads' else ((threads, v, extra) if f == '--ids' else (threads, ids, v))
binp, copies, chunks, tf, name = a[:5]
truth = dict(l.split()[:2] for l in open(tf) if l.strip())
names = [l[1:].split()[0] for l in open(copies) if l.startswith('>')]
BIG = 10
def score(groups):
    pur = sum(collections.Counter(truth[m] for m in g).most_common(1)[0][1] for g in groups.values()); return pur / sum(len(g) for g in groups.values())
def row(label, groups, sec=''):
    sizes = sorted(len(g) for g in groups.values()); big = {k: g for k, g in groups.items() if len(g) >= BIG}
    cov = sum(len(g) for g in big.values()) / len(names)
    print('%-28s %6d %6d %7d | %5.1f%% %4d %s | %.3f | %s' % (label, len(groups), sizes[len(sizes) // 2], sum(s == 1 for s in sizes), 100 * cov, len(big), ('%.3f' % score(big)) if big else '  nan', score(groups), sec), flush=True)
print('%-28s %6s %6s %7s | %s | %s | %s' % ('method', 'groups', 'med.sz', 'single', 'copies covered, n groups, purity %s (groups >= %d)' % (name, BIG), 'purity all', 'wall s'))
row('one family-wide consensus', {'all': names})
g = collections.defaultdict(list)
for l in open(chunks):
    f = l.split()
    if len(f) >= 2: g[f[1]].append(f[0])
row('SubFam chunks (%d)' % len(g), g)
with tempfile.TemporaryDirectory() as t:
    for i in ids.split(','):
        out = os.path.join(t, 'c%s.clstr' % i); t0 = time.time()
        r = subprocess.run([binp, '--id', i, '--threads', str(threads), '--output', out] + extra.split() + [copies], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sec = '%.0f' % (time.time() - t0)
        if r.returncode or not os.path.exists(out): print('%-28s failed (exit %d)' % ('meshclust2 --id ' + i, r.returncode)); continue
        g = collections.defaultdict(list); cur = None
        for l in open(out):
            if l.startswith('>Cluster'): cur = l.split()[1]
            elif cur is not None and '>' in l: g[cur].append(l.split('>', 1)[1].split('...')[0].split()[0])
        missing = len(names) - sum(len(x) for x in g.values())
        row('meshclust2 --id %s' % i + (' (%d copies missing)' % missing if missing else ''), g, sec)
