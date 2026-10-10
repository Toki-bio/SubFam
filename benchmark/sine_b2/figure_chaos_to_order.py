#!/usr/bin/env python3
"""figure_chaos_to_order.py SUB.fa SUB.tsv CHUNKS.tsv CONS.fa OUT.png [--threads T]

Sketch of the 'chaos to order' figure. SUB.fa: a hand-picked set of real B2 copies; SUB.tsv: copy<TAB>label
(B2a-e best hit); CHUNKS.tsv, CONS.fa: SubFam output on SUB.fa. All copies and all chunk consensuses are
aligned together (MAFFT --auto) so the three panels share their columns:
  a) the copies in random (input) order, b) the same copies in SubFam's k-mer order, c) the chunk consensuses.
Colour strip: the label of each copy (a, b) / the majority label of the chunk (c)."""
import collections, subprocess, sys, tempfile, os
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

args = sys.argv[1:]; threads = 4
if '--threads' in args: i = args.index('--threads'); threads = int(args[i + 1]); del args[i:i + 2]
sub, tsv, chunks, cons, out = args
def rd(p):
    d, n = collections.OrderedDict(), None
    for l in open(p):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l.strip())
    return collections.OrderedDict((k, ''.join(v).upper()) for k, v in d.items())
copies, cs = rd(sub), rd(cons)
lab = dict(l.split()[:2] for l in open(tsv))
order = [l.split()[0] for l in open(chunks) if l.strip()]        # SubFam order
chunk_of = {l.split()[0]: l.split()[1] for l in open(chunks) if l.strip()}
rows = collections.OrderedDict(copies); rows.update({'CHUNK_' + k: v for k, v in cs.items()})
with tempfile.TemporaryDirectory() as t:
    with open(t + '/in.fa', 'w') as o:
        for k, v in rows.items(): o.write('>%s\n%s\n' % (k, v))
    aln = subprocess.run(['mafft', '--auto', '--thread', str(threads), '--quiet', t + '/in.fa'], capture_output=True, text=True, check=True).stdout
A = rd_aln = collections.OrderedDict(); n = None
for l in aln.splitlines():
    if l.startswith('>'): n = l[1:].split()[0]; A[n] = []
    elif n: A[n].append(l.strip())
A = {k: ''.join(v).upper() for k, v in A.items()}
L = len(next(iter(A.values())))
code = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
def mat(names): return np.array([[code.get(c, 0) for c in A[k]] for k in names])
# drop columns that are gaps in nearly everything (alignment artefacts of the few long insertions stay visible)
keep = [j for j in range(L) if sum(A[k][j] != '-' for k in copies) >= 0.25 * len(copies)]
# colour = agreement with the column majority of ALL copies: 0 gap, 1 matches the majority, 2 differs from it
maj_col = {j: collections.Counter(A[k][j] for k in copies if A[k][j] != '-').most_common(1)[0][0] for j in keep}
def m(names): return np.array([[0 if A[k][j] == '-' else (1 if A[k][j] == maj_col[j] else 2) for j in keep] for k in names])
rnd = list(copies); ordered = order
byc = collections.defaultdict(list)
for k in ordered: byc[chunk_of[k]].append(lab[k])
chs = sorted(cs); maj = {c: collections.Counter(byc[c]).most_common(1)[0] for c in chs}
labs = ['B2a', 'B2b', 'B2c', 'B2d', 'B2e']
col = {'B2a': '#1b7837', 'B2b': '#7fbf7b', 'B2c': '#e08214', 'B2d': '#762a83', 'B2e': '#9e9e9e'}
base = ListedColormap(['#ffffff', '#dfe7ee', '#c0392b'])  # gap, matches the column majority, differs from it
fig = plt.figure(figsize=(11, 12))
gs = fig.add_gridspec(3, 2, width_ratios=[1, 28], height_ratios=[len(rnd), len(rnd), len(chs) * 9], hspace=0.12, wspace=0.01)
def panel(r, names, title, labelof):
    ax = fig.add_subplot(gs[r, 1]); ax.imshow(m(names), aspect='auto', cmap=base, vmin=0, vmax=2, interpolation='nearest')
    ax.set_title(title, loc='left', fontsize=10); ax.set_yticks([]); ax.set_xlabel('alignment column' if r == 2 else '')
    if r < 2: ax.set_xticks([])
    sb = fig.add_subplot(gs[r, 0]); sb.imshow(np.array([[labs.index(labelof(k)) if labelof(k) in labs else 4] for k in names]), aspect='auto',
        cmap=ListedColormap([col[x] for x in labs]), vmin=0, vmax=4, interpolation='nearest'); sb.axis('off')
    return ax
panel(0, rnd, 'a  %d real B2 copies in input (random) order: B2a 120, B2b 30, B2c 120, B2d 120, picked by the 2021 best-hit label' % len(rnd), lambda k: lab[k])
panel(1, ordered, 'b  the same copies, in SubFam k-mer tree order (nothing else changed)', lambda k: lab[k])
cn = ['CHUNK_' + c for c in chs]
ax = panel(2, cn, 'c  %d chunk consensuses (N = 20), the alignment SubFam passes downstream' % len(chs), lambda k: maj[k[6:]][0])
ax.yaxis.tick_right(); ax.set_yticks(range(len(chs))); ax.set_yticklabels(['%s  %s p=%.2f' % (c, maj[c][0], maj[c][1] / len(byc[c])) for c in chs], fontsize=7)
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col[x]) for x in labs] + [plt.Rectangle((0, 0), 1, 1, color=c) for c in base.colors[1:]],
           labels=labs + ['= majority base of the column', 'differs from it (white = gap)'], ncol=7, loc='lower center', fontsize=8, frameon=False)
fig.savefig(out, dpi=170, bbox_inches='tight'); print('wrote', out, len(rnd), 'copies', len(chs), 'chunks', L, 'columns,', len(keep), 'shown', file=sys.stderr)
