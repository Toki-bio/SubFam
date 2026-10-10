#!/usr/bin/env python3
"""check_rows.py EXAMPLE.aln.fasta  -  identity of every SubFam row to the head of its group.

Rows are grouped as in the example files: a head row (TRUE_... or PUBLISHED_...) followed by the
SubFam rows of that group. Identity = matching columns / columns where both have a residue, over
the span where both have residues (end gaps excluded). Prints one line per row and a summary
(median, min, and how many rows are closer to another group's head).
"""
import statistics, sys


def rd(p):
    out, n = [], None
    for l in open(p):
        if l.startswith('>'): n = l[1:].strip(); out.append([n, []])
        elif n is not None: out[-1][1].append(l.strip())
    return [(n, ''.join(s).upper()) for n, s in out]


def ident(a, b):
    m = t = 0
    for x, y in zip(a, b):
        if x != '-' and y != '-':
            t += 1; m += x == y
    return m / t if t else 0.0


rows = [r for r in rd(sys.argv[1]) if not r[0].startswith('REFERENCE_')]  # REFERENCE_ rows belong to no group
is_head = [n.startswith(('TRUE_', 'PUBLISHED_', 'CONS_')) for n, _ in rows]
# a group's head is the first head row of a block of head rows (TRUE_ followed by PUBLISHED_ with --both)
heads = [(i, n, s) for i, (n, s) in enumerate(rows) if is_head[i] and not (i and is_head[i - 1])]
ids, closer = [], 0
for i, (n, s) in enumerate(rows):
    if is_head[i]: continue
    own = max((h for h in heads if h[0] < i), key=lambda h: h[0], default=None)
    if own is None: continue
    d = ident(s, own[2]); ids.append(d)
    best = max(heads, key=lambda h: ident(s, h[2]))
    flag = '' if best[0] == own[0] else '  closer to ' + best[1]
    if flag: closer += 1
    print(f'{n}\t{own[1]}\t{d:.3f}{flag}')
if ids:
    print(f'# {len(ids)} rows: identity to own head median {statistics.median(ids):.3f} '
          f'min {min(ids):.3f} max {max(ids):.3f}; {closer} rows closer to another head', file=sys.stderr)
