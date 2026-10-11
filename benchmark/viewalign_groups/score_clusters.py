#!/usr/bin/env python3
"""score_clusters.py CLUSTERS.json CHUNKS.tsv NAME1 TRUTH1.tsv [NAME2 TRUTH2.tsv ...]
Prints the cluster hierarchy of cluster_rows.js; for every cluster: size, number of exclusive (perfect) columns, whether its rows are
consecutive chunks in SubFam's order (the rows are named <chunk> or <chunk>|..., chunk ids sort in order), and the majority label of
the rows' chunks for every labelling, with the purity (fraction of the cluster's COPIES carrying that label, i.e. chunk members counted)."""
import collections, json, sys
js, chunks = sys.argv[1:3]; lab = []
for i in range(3, len(sys.argv) - 1, 2): lab.append((sys.argv[i], dict(l.split()[:2] for l in open(sys.argv[i + 1]) if l.strip())))
members = collections.defaultdict(list)
for l in open(chunks):
    f = l.split()
    if len(f) >= 2: members[f[1]].append(f[0])
order = {c: i for i, c in enumerate(sorted(members))}
d = json.load(open(js))
def show(c, ind):
    rows = [r.split('|')[0] for r in c['rows']]; pos = sorted(order[r] for r in rows)
    runs = 1 + sum(b - a > 1 for a, b in zip(pos, pos[1:]))
    s = '%s%d rows, %d exclusive columns, %d run(s) in chunk order (%s..%s)' % (ind, c['size'], len(c['perfect']), runs, rows and min(rows, key=order.get), rows and max(rows, key=order.get))
    for name, t in lab:
        cnt = collections.Counter(t[m] for r in rows for m in members[r] if m in t); n = sum(cnt.values())
        s += ' | %s %s' % (name, ', '.join('%s %d' % (k, v) for k, v in cnt.most_common(3)))
        s += ' (%.2f)' % (cnt.most_common(1)[0][1] / n)
    print(s)
    for ch in c.get('children', []): show(ch, ind + '    ')
    if c.get('childUnassigned'): print(ind + '    unassigned inside: %d rows' % len(c['childUnassigned']))
for c in d['clusters']: show(c, '')
if d['unassigned']: print('unassigned at top level: %d rows' % len(d['unassigned']))
