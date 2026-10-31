#!/usr/bin/env python3
"""figure_grouped_alignment.py GROUPED.aln.fasta OUT.png|OUT.pdf [--title T] [--min-occ 0.25] [--share] [--no-refs]
Figure 1: a grouped example alignment (examples/*.aln.fasta, rows = CLASSCONS_/PUBLISHED_ heads followed by SubFam rows) drawn as a
raster: one row per sequence, columns present in >= min-occ of the SubFam rows, colour = agreement with the column majority of the SubFam
rows (grey = majority base, red = other base, white = gap); head rows in a darker frame; a left strip colours the group; purity of each
SubFam row as text on the right only with --share. REFERENCE_ rows (ungrouped) are drawn at the top unless --no-refs. Needs matplotlib."""
import re, sys, collections
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
a = sys.argv[1:]; title = ''; minocc = 0.25
share = '--share' in a; norefs = '--no-refs' in a; a = [x for x in a if x not in ('--share', '--no-refs')]
if '--title' in a: i = a.index('--title'); title = a[i + 1]; del a[i:i + 2]
if '--min-occ' in a: i = a.index('--min-occ'); minocc = float(a[i + 1]); del a[i:i + 2]
inp, out = a[:2]
rows = []; n = None
for l in open(inp):
    if l.startswith('>'): n = l[1:].strip(); rows.append([n, []])
    elif n is not None: rows[-1][1].append(l.strip())
rows = [(k, ''.join(v).upper()) for k, v in rows]
if norefs: rows = [(k, v) for k, v in rows if not k.startswith('REFERENCE_')]
sub = [s for k, s in rows if '|' in k]; L = len(rows[0][1])
keep = [j for j in range(L) if sum(s[j] != '-' for s in sub) >= minocc * len(sub)]
maj = {j: collections.Counter(s[j] for s in sub if s[j] != '-').most_common(1)[0][0] for j in keep}
M = np.array([[0 if s[j] == '-' else (1 if s[j] == maj[j] else 2) for j in keep] for _, s in rows])
# groups: a head block (CLASSCONS_/PUBLISHED_) starts a group; REFERENCE_ rows are group -1
grp = []; g = -1; prev_head = False
for k, _ in rows:
    if k.startswith('REFERENCE_'): grp.append(-1); prev_head = False; continue
    head = k.startswith(('CLASSCONS_', 'PUBLISHED_', 'TRUE_', 'CONS_'))
    if head and not prev_head: g += 1
    grp.append(g); prev_head = head
names = []
for k, _ in rows:
    if k.startswith('REFERENCE_'): names.append(k.replace('REFERENCE_', 'ref ')) 
    elif '|' in k: f = k.split('|'); names.append('%s  %s' % (f[0], f[2] if len(f) > 2 else '') if share else f[0])
    else: names.append(re.sub(r'_consensus_of_.*', '', k).replace('CLASSCONS_', 'class consensus ').replace('PUBLISHED_', 'published '))
if g < 0:                      # no head rows: group = the class named in the row (<chunk>|<class>|p<purity>), in order of first appearance
    classes = []
    for k, _ in rows:
        c = k.split('|')[1] if '|' in k else '?'
        if c not in classes: classes.append(c)
    grp = [classes.index(k.split('|')[1]) if '|' in k else -1 for k, _ in rows]; g = len(classes) - 1
ngr = g + 1; pal = plt.get_cmap('tab10' if ngr <= 10 else 'tab20')
fig_h = max(4, min(14, 0.16 * len(rows) + 1.2))
fig = plt.figure(figsize=(12, fig_h)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 40], wspace=0.005)
ax = fig.add_subplot(gs[0, 1]); ax.imshow(M, aspect='auto', cmap=ListedColormap(['#ffffff', '#dfe7ee', '#c0392b']), vmin=0, vmax=2, interpolation='nearest')
ax.yaxis.tick_right()
if len(rows) <= 130: ax.set_yticks(range(len(rows))); ax.set_yticklabels(names, fontsize=6, family='monospace')
else: ax.set_yticks([]); ax.set_ylabel('%d rows, in SubFam order' % len(rows))
for i, (k, _) in enumerate(rows):
    if not ('|' in k or k.startswith('REFERENCE_')): ax.axhspan(i - 0.5, i + 0.5, facecolor='none', edgecolor='#333', lw=0.8)
ax.set_xlabel('alignment column (columns present in >= %d %% of the SubFam rows)' % round(100 * minocc)); ax.set_title(title, loc='left', fontsize=10)
sb = fig.add_subplot(gs[0, 0]); sb.imshow(np.array([[x if x >= 0 else ngr] for x in grp]), aspect='auto', cmap=ListedColormap([pal(i) for i in range(ngr)] + ['#bbbbbb']), vmin=0, vmax=ngr, interpolation='nearest'); sb.axis('off')
if 'classes' in dir():
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=pal(i)) for i in range(ngr)], labels=classes, ncol=min(ngr, 6), loc='lower center', fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.savefig(out, dpi=300, bbox_inches='tight'); print('wrote', out, len(rows), 'rows,', len(keep), 'of', L, 'columns', file=sys.stderr)
