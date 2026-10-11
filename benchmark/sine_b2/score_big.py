#!/usr/bin/env python3
"""score_big.py WORKDIR RUN [RUN ...]
B2 stress case (mm39 chr1-3, 30,000 copies): for every SubFam run in WORKDIR/RUN (prefix b2): rows, distinct rows, chunk purity against the
RepeatMasker class (id prefix), mean 'agree' of chunkstats, recovery of the five Dfam B2 consensuses (../sine_b2 dfam set: WORKDIR/dfam_B2.fa) and
of the 60 chr19 rows (WORKDIR/chr19_rows.fa) by ../alu_hg38/recover.py (coverage >= 0.9; identity 0.95/0.98/0.99), and the reverse view."""
import collections, os, subprocess, sys
w = sys.argv[1]; runs = sys.argv[2:]; here = os.path.dirname(os.path.abspath(__file__)); rec = os.path.join(here, '..', 'alu_hg38', 'recover.py')
sample = [l[1:].split()[0] for l in open(os.path.join(w, 's1.fa')) if l.startswith('>')]; S = len(sample)
truth = {n: n.split('_chr')[0] for n in sample}
def rd(p):
    d = collections.OrderedDict(); n = None
    for l in open(p):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}
for run in runs:
    d = os.path.join(w, run); cons = rd(os.path.join(d, 'b2.cons.fasta'))
    m = collections.defaultdict(list)
    for l in open(os.path.join(d, 'b2.chunks.tsv')):
        f = l.split(); m[f[1]].append(f[0])
    pur = sum(collections.Counter(truth[x] for x in g).most_common(1)[0][1] for g in m.values()) / S
    per = collections.Counter(collections.Counter(truth[x] for x in g).most_common(1)[0][0] for g in m.values())
    st = [l.split('\t') for l in open(os.path.join(d, 'b2.chunkstats.tsv')) if not l.startswith('chunk')]
    print('== %s: %d rows, %d distinct, chunk purity vs rmsk class %.3f, mean agree %.3f; rows per majority class %s' % (run, len(cons), len(set(cons.values())), pur, sum(float(r[5]) for r in st) / len(st), dict(per)), flush=True)
    for ref in ('dfam_B2.fa', 'chr19_rows.fa'):
        out = os.path.join(d, 'rec_' + ref.split('.')[0])
        if not os.path.exists(out + '.ref_best.tsv'):
            subprocess.run(['python3', rec, os.path.join(w, ref), os.path.join(d, 'b2.cons.fasta'), out], check=True, stdout=subprocess.DEVNULL)
        rows = [l.rstrip('\n').split('\t') for l in open(out + '.ref_best.tsv')][1:]; v = [float(r[3]) for r in rows]
        print('   %s (%d references): recovered at 0.95/0.98/0.99: %d/%d/%d; median best identity %.3f; min %.3f' % (ref, len(v), sum(x >= .95 for x in v), sum(x >= .98 for x in v), sum(x >= .99 for x in v), sorted(v)[len(v) // 2], min(v)))
        if ref == 'dfam_B2.fa':
            for r in rows: print('      %-22s best row %-8s identity %.3f coverage %s' % (r[0], r[2], float(r[3]), r[4]))
        cb = [l.rstrip('\n').split('\t') for l in open(out + '.cons_best.tsv')][1:]
        near = sum(float(r[2]) >= 0.95 for r in cb); print('      rows with a reference at >= 0.95: %d of %d' % (near, len(cb)), flush=True)
