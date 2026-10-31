#!/usr/bin/env python3
"""figure_chaos_to_order.py SUB.fa SUB.tsv CHUNKS.tsv CONS.fa OUT.png [--threads T]

Figure 2, 'chaos to order' (OUT.png at 300 dpi and OUT.pdf, 17 cm wide). SUB.fa: a stratified random sample of real B2 copies (chaos_subset.py); SUB.tsv: copy<TAB>label
(B2a-e best hit); CHUNKS.tsv, CONS.fa: SubFam output on SUB.fa. All copies and all chunk consensuses are
aligned together (MAFFT --auto) so the three panels share their columns:
  a) the copies in random (input) order, b) the same copies in SubFam's k-mer order, c) the chunk consensuses.
Colour strip: the label of each copy (a, b) / the majority label of the chunk (c)."""
import collections, subprocess, sys, tempfile, os
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 7, 'pdf.fonttype': 42})
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
    # one thread: the multithreaded alignment differed between runs (1,106 and 1,100 columns)
    aln = subprocess.run(['mafft', '--auto', '--thread', '1', '--quiet', t + '/in.fa'], capture_output=True, text=True, check=True).stdout
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
fig = plt.figure(figsize=(6.7, 8.4))
gs = fig.add_gridspec(3, 2, width_ratios=[1, 28], height_ratios=[len(rnd), len(rnd), len(chs) * 12], hspace=0.12, wspace=0.01)
def panel(r, names, title, labelof):
    ax = fig.add_subplot(gs[r, 1]); ax.imshow(m(names), aspect='auto', cmap=base, vmin=0, vmax=2, interpolation='nearest')
    ax.set_title(title[3:], loc='left', fontsize=7.5); ax.text(-0.075, 1.02, title[0], transform=ax.transAxes, fontweight='bold', fontsize=9, va='bottom'); ax.set_yticks([]); ax.set_xlabel('alignment column' if r == 2 else '')
    if r < 2: ax.set_xticks([])
    sb = fig.add_subplot(gs[r, 0]); sb.imshow(np.array([[labs.index(labelof(k)) if labelof(k) in labs else 4] for k in names]), aspect='auto',
        cmap=ListedColormap([col[x] for x in labs]), vmin=0, vmax=4, interpolation='nearest'); sb.axis('off')
    return ax
panel(0, rnd, 'a  %d B2 copies, input order' % len(rnd), lambda k: lab[k])
panel(1, ordered, 'b  the same copies in SubFam order', lambda k: lab[k])
cn = ['CHUNK_' + c for c in chs]
ax = panel(2, cn, 'c  their %d chunk consensuses (N = 20)' % len(chs), lambda k: maj[k[6:]][0])
ax.yaxis.tick_right(); ax.set_yticks(range(len(chs))); ax.set_yticklabels(['%d %s' % (i + 1, maj[c][0]) for i, c in enumerate(chs)], fontsize=5)
used = [x for x in labs if any(lab[k] == x for k in rnd)]   # only the groups present in the sample
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col[x]) for x in used] + [plt.Rectangle((0, 0), 1, 1, color=c, ec='#999999', lw=0.3) for c in base.colors],
           labels=used + ['gap', 'majority base of the column', 'other base'], ncol=7, loc='lower center', fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.02))
fig.savefig(out, dpi=300, bbox_inches='tight'); fig.savefig(os.path.splitext(out)[0] + '.pdf', bbox_inches='tight'); print('wrote', out, len(rnd), 'copies', len(chs), 'chunks', L, 'columns,', len(keep), 'shown', file=sys.stderr)
