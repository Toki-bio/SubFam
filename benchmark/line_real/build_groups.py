#!/usr/bin/env python3
"""Real L1 copies: SubFam chunk consensuses -> peel groups (SINE-discriminator peel_features.py) -> one consensus per group
built from the ORIGINAL copies of the group's chunks (not from the chunk consensuses) -> one alignment, grouped.

usage: build_groups.py COPIES.fa SUBFAM_DIR PREFIX PEEL.json OUT.aln.fasta [--sample N] [--threads T] [--tmp DIR]
  SUBFAM_DIR/PREFIX.cons.fasta, PREFIX.chunks.tsv   PEEL.json = peel_features.json (key "peeled": list of {members:[chunk,...]})
Rows: CONS_G<k>_<n>copies (group head), then <chunk>|G<k>|<rmsk majority>|p<purity of that chunk>; chunks not peeled: G0 (residue).
A group consensus uses at most N copies (random, seed 1), because MAFFT on hundreds of multi-kb copies is slow; the sample size is in the row name."""
import sys, os, json, random, collections, argparse, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'examples'))
from build_examples import rd, mafft, plurality

ap = argparse.ArgumentParser()
ap.add_argument('copies'); ap.add_argument('sfdir'); ap.add_argument('prefix'); ap.add_argument('peel'); ap.add_argument('out')
ap.add_argument('--sample', type=int, default=60); ap.add_argument('--threads', type=int, default=4); ap.add_argument('--tmp', default='/tmp/bg')
a = ap.parse_args(); os.makedirs(a.tmp, exist_ok=True)
copies = rd(a.copies); cons = rd(os.path.join(a.sfdir, a.prefix + '.cons.fasta'))
chunk_copies = collections.defaultdict(list)
for l in open(os.path.join(a.sfdir, a.prefix + '.chunks.tsv')):
    f = l.split()
    if len(f) >= 2: chunk_copies[f[1]].append(f[0])
lab = lambda c: c.split('_')[0]
peel = json.load(open(a.peel))['peeled']
group_of = {}; order = []
for k, g in enumerate(peel, 1):
    for ch in g['members']: group_of[ch] = k
    order.append(k)
if any(ch not in group_of for ch in cons): order.append(0)
rows = collections.OrderedDict(); rnd = random.Random(1); report = []
for k in order:
    chunks = [ch for ch in cons if group_of.get(ch, 0) == k]
    cps = [c for ch in chunks for c in chunk_copies[ch]]
    labs = collections.Counter(lab(c) for c in cps)
    if k and len(cps) >= 3:
        pick = cps if len(cps) <= a.sample else rnd.sample(cps, a.sample)
        al = mafft({c: copies[c] for c in pick}, a.threads, a.tmp, 'g%d' % k, ('--retree', '2', '--maxiterate', '0'))
        rows['CONS_G%d_%dof%dcopies' % (k, len(pick), len(cps))] = plurality(al)
    for ch in chunks:
        cl = collections.Counter(lab(c) for c in chunk_copies[ch]); m, n = cl.most_common(1)[0]
        rows['%s|G%d|%s|p%.2f' % (ch, k, m, n / len(chunk_copies[ch]))] = cons[ch]
    report.append((k, len(chunks), len(cps), dict(labs.most_common(4))))
    print('G%d: %d chunks, %d copies, rmsk %s' % report[-1], file=sys.stderr)
al = mafft(rows, a.threads, a.tmp, 'final', ('--auto',))
with open(a.out, 'w') as o:
    for k in rows: o.write('>%s\n%s\n' % (k, al[k]))
print('wrote', a.out, len(rows), 'rows', file=sys.stderr)
