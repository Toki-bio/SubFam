#!/usr/bin/env python3
"""chaos_subset.py COPIES.fa LABELS.tsv OUT.fa OUT.tsv
The hand-picked set of the chaos-to-order figure: a random (seed 1) 120 B2a, 30 B2b, 120 B2c and 120 B2d copies
(labels = best hit to the Vassetzky et al. 2021 consensuses), written in shuffled order (no B2e: its copies
are not separable from B2d at the copy level, see README)."""
import collections, random, sys
fa, tsv, out, outt = sys.argv[1:5]
lab = dict(l.split()[:2] for l in open(tsv)); seqs = collections.OrderedDict(); n = None
for l in open(fa):
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []
    else: seqs[n].append(l.strip())
by = collections.defaultdict(list)
for k, v in lab.items(): by[v].append(k)
r = random.Random(1); pick = []
for sf, m in (('B2a', 120), ('B2b', 30), ('B2c', 120), ('B2d', 120)): pick += r.sample(sorted(by[sf]), m)
r.shuffle(pick)
with open(out, 'w') as o:
    for k in pick: o.write('>%s\n%s\n' % (k, ''.join(seqs[k])))
open(outt, 'w').write(''.join('%s\t%s\n' % (k, lab[k]) for k in pick))
