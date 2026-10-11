#!/usr/bin/env python3
"""recompute_3end.py: the numbers of README sections 2-3 from the committed raw files.
  hits3p.tsv.gz        best hits of the last 900 bp of every copy to the Dfam L1PA2-8 3'-end consensuses (relabel_3end.sh)
  sf_last900/*.chunks.tsv   SubFam -n 20 on the last 900 bp (l1_3p.fa); ../../benchmark/baseline_vsearch l1_sweep.txt holds the whole-copy run
Prints: relabel agreement with rmsk; chunk purity and chunks >= 0.9 for the 3'-end run and, if given as argv[1], a whole-copy chunks.tsv."""
import collections, glob, gzip, sys
hits = collections.defaultdict(list)
for l in gzip.open('hits3p.tsv.gz', 'rt'):
    q, t, i, qc, tc = l.split(); hits[q].append((float(i), t.split('_')[0] if t.startswith('L1PA') else t.split('|')[0]))
lab = {}
for q, h in hits.items():
    best = max(x[0] for x in h); tops = {t for i, t in h if i == best}; lab[q] = tops
same = uniq = within1 = tie = 0; n = len(lab)
order = ['L1PA%d' % k for k in range(2, 9)]
for q, tops in lab.items():
    r = q.split('_')[0]; tie += len(tops) > 1
    if r in tops: same += 1
    if tops == {r}: uniq += 1
    if any(abs(order.index(r) - order.index(t)) <= 1 for t in tops if t in order): within1 += 1
print('copies with a 3\'-end hit: %d; unique best hit equal to rmsk: %d (%.1f %%); rmsk among the best hits (ties counted): %d (%.1f %%); within one class: %d (%.1f %%); ties: %d' % (n, uniq, 100 * uniq / n, same, 100 * same / n, within1, 100 * within1 / n, tie))
def purity(chunks):
    m = collections.defaultdict(list)
    for l in open(chunks):
        f = l.split(); m[f[1]].append(f[0].split('_')[0])
    p = [collections.Counter(v).most_common(1)[0][1] / len(v) for v in m.values()]
    return len(p), sum(collections.Counter(v).most_common(1)[0][1] for v in m.values()) / sum(len(v) for v in m.values()), sum(x >= 0.9 for x in p)
for name, f in [('last 900 bp', glob.glob('sf_last900/*.chunks.tsv')[0])] + ([('whole copies', sys.argv[1])] if len(sys.argv) > 1 else []):
    k, pur, hi = purity(f); print('%s: %d chunks, purity %.3f, chunks with purity >= 0.9: %d' % (name, k, pur, hi))
