# Compares COSEG assignments and SubFam chunks against the lineage labels (loci_labels.tsv), ambiguous labels removed.
# usage: compare.py loci_labels.tsv names.txt NAME  [coseg_dir:label ...] [subfam_chunks.tsv:label ...]
# Metrics: groups, purity (majority label share), homogeneity, completeness, V-measure.
import sys, collections, math
labf, namesf = sys.argv[1:3]
lab = {l.split('\t')[0]: l.split('\t')[6] for l in open(labf) if not l.startswith('acc')}
names = [l.strip().split('|')[0] for l in open(namesf)]
def H(c):
    n = sum(c.values()); return -sum(v / n * math.log(v / n) for v in c.values())
def vm(cl, tr):
    n = len(cl); ct = collections.Counter(zip(cl, tr)); bycl = collections.defaultdict(collections.Counter); bytr = collections.defaultdict(collections.Counter)
    for (k, t), v in ct.items(): bycl[k][t] += v; bytr[t][k] += v
    Hc, Hk = H(collections.Counter(tr)), H(collections.Counter(cl))
    hom = 1 - sum(sum(c.values()) / n * H(c) for c in bycl.values()) / Hc if Hc else 1
    comp = 1 - sum(sum(c.values()) / n * H(c) for c in bytr.values()) / Hk if Hk else 1
    return hom, comp, 2 * hom * comp / (hom + comp)
def pur(cl, tr):
    c = collections.defaultdict(collections.Counter)
    for k, t in zip(cl, tr): c[k][t] += 1
    return sum(max(v.values()) for v in c.values()) / len(cl)
idx = [i for i, n in enumerate(names) if lab[n] != 'ambig']
tr = [lab[names[i]] for i in idx]
print('%-22s %6s %6s %6s %6s %6s' % ('method', 'groups', 'purity', 'homog', 'compl', 'V'))
for spec in sys.argv[3:]:
    path, label = spec.rsplit(':', 1)
    if path.endswith('.assign'):
        asg = [l.split()[-1] for l in open(path)]                      # one line per copy, in input order
    else:
        ch = {l.split('\t')[0].split('|')[0]: l.split('\t')[1] for l in open(path)}
        asg = [ch[n] for n in names]
    cl = [asg[i] for i in idx]; h, c, v = vm(cl, tr)
    print('%-22s %6d %6.3f %6.3f %6.3f %6.3f' % (label, len(set(cl)), pur(cl, tr), h, c, v))
