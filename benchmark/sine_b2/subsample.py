#!/usr/bin/env python3
"""subsample.py IN.fa N SEED OUT.fa TRUTH.tsv
Random N records of IN.fa (python random, the given seed), in input order; TRUTH.tsv maps each
id to its RepeatMasker name, i.e. the id without the trailing _chrN_start_end."""
import random, re, sys

inp, n, seed, out, truth = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
names, seqs = [], []
for l in open(inp):
    if l.startswith('>'): names.append(l[1:].strip()); seqs.append([])
    else: seqs[-1].append(l.strip())
idx = list(range(len(names))); random.Random(seed).shuffle(idx); idx = sorted(idx[:n])
with open(out, 'w') as o, open(truth, 'w') as t:
    for i in idx:
        o.write('>%s\n%s\n' % (names[i], ''.join(seqs[i])))
        t.write('%s\t%s\n' % (names[i], re.sub(r'_chr[^_]+_\d+_\d+$', '', names[i])))
print(len(idx), 'of', len(names), file=sys.stderr)
