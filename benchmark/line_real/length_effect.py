#!/usr/bin/env python3
"""length_effect.py L1.fa CHUNKS.tsv    (label = the part of the id before the first '_')
Is SubFam's k-mer order driven by copy length rather than by sequence differences? (the weighted Jaccard distance of
k-mer counts is large between a sequence and its own prefix). Prints: (1) the distance of the longest copy to its own
prefixes; (2) normalised mutual information of the chunk with the length quintile and with the label; (3) purity of
chunks cut by length alone; within-chunk length spread against random chunks; (4) label purity of the chunks inside each length quintile."""
import collections, math, random, statistics as st, sys
fa, chunks = sys.argv[1:3]
seqs, n = {}, None
for l in open(fa):
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []
    elif n: seqs[n].append(l.strip())
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}; ids = list(seqs)
lab = {k: k.split('_')[0] for k in ids}; L = {k: len(v) for k, v in seqs.items()}
ch = {l.split()[0]: l.split()[1] for l in open(chunks) if l.strip()}
kc = lambda s, k=6: collections.Counter(s[i:i + k] for i in range(len(s) - k + 1))
def wj(a, b): ks = set(a) | set(b); return sum(min(a[x], b[x]) for x in ks) / sum(max(a[x], b[x]) for x in ks)
k0 = max(ids, key=L.get); print('(1) longest copy %d bp; weighted-Jaccard distance to its own prefix of:' % L[k0])
for m in (3000, 2000, 1000, 500): print('    %4d bp: %.3f' % (m, 1 - wj(kc(seqs[k0]), kc(seqs[k0][:m]))))
def H(a): return -sum(v / len(a) * math.log(v / len(a)) for v in collections.Counter(a).values())
def nmi(a, b):
    n = len(a); ca, cb, cab = collections.Counter(a), collections.Counter(b), collections.Counter(zip(a, b))
    return sum(v / n * math.log((v / n) / ((ca[x] / n) * (cb[y] / n))) for (x, y), v in cab.items()) / math.sqrt(H(a) * H(b))
q = sorted(L.values()); cut = [q[len(q) * i // 5] for i in (1, 2, 3, 4)]
lq = [sum(L[k] > c for c in cut) for k in ids]; labs = [lab[k] for k in ids]; chs = [ch[k] for k in ids]
print('(2) NMI(chunk; length quintile) %.3f   NMI(chunk; label) %.3f   NMI(label; length quintile) %.3f' % (nmi(chs, lq), nmi(chs, labs), nmi(labs, lq)))
purity = lambda gs: sum(collections.Counter(lab[m] for m in g).most_common(1)[0][1] for g in gs) / len(ids)
by = sorted(ids, key=L.get); gl = [by[i:i + 20] for i in range(0, len(by), 20)]
gs = [[k for k in ids if ch[k] == c] for c in sorted(set(chs))]
sh = ids[:]; random.Random(1).shuffle(sh); gr = [sh[i:i + 20] for i in range(0, len(sh), 20)]
cv = lambda G: st.median(st.pstdev([L[k] for k in g]) / st.mean([L[k] for k in g]) for g in G)
print('(3) purity: chunks cut by length alone %.3f, SubFam chunks %.3f; median within-chunk length CV: SubFam %.2f, random chunks %.2f' % (purity(gl), purity(gs), cv(gs), cv(gr)))
print('(4) label purity of SubFam chunks inside each length quintile (chunks with at least 5 copies of the quintile):')
for b in range(5):
    byc = collections.defaultdict(list)
    for k, x in zip(ids, lq):
        if x == b: byc[ch[k]].append(lab[k])
    big = [v for v in byc.values() if len(v) >= 5]; tot = sum(len(v) for v in big)
    print('    quintile %d (%d-%d bp): %d copies in %d such chunks, purity %.3f' % (b + 1, min(L[k] for k, x in zip(ids, lq) if x == b), max(L[k] for k, x in zip(ids, lq) if x == b), tot, len(big), sum(collections.Counter(v).most_common(1)[0][1] for v in big) / tot))
