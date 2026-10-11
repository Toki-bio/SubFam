#!/usr/bin/env python3
"""region_topology.py L1.fa OUT_PREFIX [--threads T]

Does the k-mer tree of the L1 copies depend on which part of the sequence it is built from?
Same copies, same distance (weighted Jaccard of 6-mer counts, UPGMA), different regions. The copies are 5'-truncated and
end at the 3' end, so the last X bp of every copy that is at least X long is one homologous window (up to indels).
  A  all copies:            whole sequence  vs  last 1000 bp (every copy has it; the shortest copy is 1,000 bp)
  B  copies >= 2000 bp:     whole sequence  vs  last 2000 bp (about a third of a full-length L1PA copy)
  C  control:               whole sequence k=6 vs k=5 on the same copies (how much a harmless parameter change moves the tree)
Compared per pair of trees (same leaves): normalised Robinson-Foulds distance of the rooted clades, cophenetic correlation,
and per tree: purity of consecutive chunks of 20 in SubFam's leaf order against the RepeatMasker class (label = id up to the first '_'),
adjusted Rand index of the tree cut into 7 clusters, NMI of the chunk with the length quintile of the whole copy.
The leaf order is SubFam's own (upgma_order of SubFam.sh, closest-ends flipping); clades, RF and cuts use scipy average linkage on the same matrix."""
import collections, math, os, re, subprocess, sys
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, cophenet
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, '..', '..', 'SubFam.sh')).read()
code = re.search(r"KMER_ORDER_PY=\$\(cat <<'PY'\n(.*?)\nPY\n\)", src, re.S).group(1)
ns = {'__name__': 'subfam_kmer'}; exec(compile(code, 'SubFam.sh:KMER_ORDER_PY', 'exec'), ns)

fa, prefix = sys.argv[1:3]
names, seqs = ns['read_fasta'](fa)
lab = [n.split('_')[0] for n in names]; L = np.array([len(s) for s in seqs])
classes = sorted(set(lab))

def tree(idx, region, k=6):
    """distance matrix -> (SubFam leaf order, scipy linkage, condensed distances) for copies idx cut to region (None = whole, X = last X bp)"""
    ss = [seqs[i][-region:] if region else seqs[i] for i in idx]
    d, profiles, total = ns['distances'](ss, k, False)
    dd = d.astype(np.float64); np.fill_diagonal(dd, 0); cond = squareform((dd + dd.T) / 2, checks=False)
    order = ns['upgma_order'](d, lambda a, b: ns['pair_distance'](profiles, total, a, b))
    return order, linkage(cond, 'average'), cond

def clades(Z, n):
    members = {i: frozenset([i]) for i in range(n)}; out = set()
    for j, (a, b, _, _) in enumerate(Z):
        c = members[int(a)] | members[int(b)]; members[n + j] = c
        if len(c) < n: out.add(c)
    return out

def rf(Z1, Z2, n):
    a, b = clades(Z1, n), clades(Z2, n); return len(a ^ b) / (len(a) + len(b))

def nmi(a, b):
    n = len(a); ca, cb, cab = collections.Counter(a), collections.Counter(b), collections.Counter(zip(a, b))
    h = lambda c: -sum(v / n * math.log(v / n) for v in c.values())
    return sum(v / n * math.log((v / n) / ((ca[x] / n) * (cb[y] / n))) for (x, y), v in cab.items()) / math.sqrt(h(ca) * h(cb))

def ari(a, b):
    from math import comb
    n = len(a); cab = collections.Counter(zip(a, b)); ca = collections.Counter(a); cb = collections.Counter(b)
    s = sum(comb(v, 2) for v in cab.values()); sa = sum(comb(v, 2) for v in ca.values()); sb = sum(comb(v, 2) for v in cb.values())
    e = sa * sb / comb(n, 2); return (s - e) / ((sa + sb) / 2 - e)

def describe(idx, order, Z, tag):
    labs = [lab[idx[i]] for i in order]; chunk = {}
    chunks = [order[i:i + 20] for i in range(0, len(order), 20)]
    for ci, g in enumerate(chunks):
        for i in g: chunk[i] = ci
    pur = sum(collections.Counter(lab[idx[i]] for i in g).most_common(1)[0][1] for g in chunks) / len(idx)
    cl = fcluster(Z, len(classes), 'maxclust')
    q = np.quantile(L, [.2, .4, .6, .8]); lq = [int(sum(L[idx[i]] > c for c in q)) for i in range(len(idx))]
    r = (tag, len(idx), len(chunks), pur, ari([lab[i] for i in idx], list(cl)), nmi([chunk[i] for i in range(len(idx))], lq))
    print('%-34s copies %4d  chunks %3d  chunk purity %.3f  ARI(7 clusters, class) %.3f  NMI(chunk, length quintile) %.3f' % r, flush=True)

allidx = list(range(len(seqs))); shared = [i for i in allidx if L[i] >= 2000]
print('copies: %d; with >= 2000 bp: %d (%.0f %%)' % (len(allidx), len(shared), 100 * len(shared) / len(allidx)), flush=True)
runs = collections.OrderedDict()
for tag, idx, region, k in (('A whole', allidx, None, 6), ('A last 1000 bp', allidx, 1000, 6), ('B whole', shared, None, 6),
                            ('B last 2000 bp', shared, 2000, 6), ('C whole, k=5', allidx, None, 5)):
    order, Z, cond = tree(idx, region, k); runs[tag] = (idx, order, Z, cond); describe(idx, order, Z, tag)
    with open('%s.%s.order.txt' % (prefix, tag.replace(' ', '_').replace(',', '').replace('=', '')), 'w') as o:
        o.write('\n'.join(names[idx[i]] for i in order) + '\n')
print()
for a, b in (('A whole', 'A last 1000 bp'), ('B whole', 'B last 2000 bp'), ('A whole', 'C whole, k=5')):
    (ia, _, Za, ca), (ib, _, Zb, cb) = runs[a], runs[b]; n = len(ia)
    print('%-18s vs %-18s  normalised RF %.3f   cophenetic correlation %.3f' % (a, b, rf(Za, Zb, n), np.corrcoef(cophenet(Za), cophenet(Zb))[0, 1]), flush=True)
