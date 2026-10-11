#!/usr/bin/env python3
"""score_runs.py WORKDIR RUN [RUN ...]   (WORKDIR holds s1.fa, s1.truth.tsv, chr1_alu_full.fa, dfam_alu.fa, price213.fa; RUN = n50, n20, n50P ...)
For every SubFam run: rows, distinct rows, chunk purity against the RepeatMasker name, mean 'agree' of chunkstats, and recovery of the
Dfam and Price consensuses (recover.py; coverage >= 0.9; identity 0.95/0.98/0.99) per eligibility stratum: expected copies of the reference
in the sample >= 50, 10-50, < 10 (Dfam: rmsk count of that name on the chromosome x sample/total; Price: his copy count x sample/his total)."""
import collections, os, subprocess, sys
w = sys.argv[1]; runs = sys.argv[2:]; here = os.path.dirname(os.path.abspath(__file__))
names = [l[1:].split()[0] for l in open(os.path.join(w, 'chr1_alu_full.fa')) if l.startswith('>')]
tot = len(names); rm = collections.Counter(n.split('_')[0] for n in names)
sample = [l[1:].split()[0] for l in open(os.path.join(w, 's1.fa')) if l.startswith('>')]; S = len(sample)
truth = {n: n.split('_')[0] for n in sample}
price = {}
for l in open(os.path.join(w, 'price213.fa')):
    if l.startswith('>'): f = l[1:].strip().split('|'); price[f[0]] = int(f[2])
ptot = sum(price.values())
def expected(ref):
    if ref.startswith('DFAM_'): return rm.get(ref[5:], 0) * S / tot
    return price.get(ref.split('|')[0], 0) * S / ptot
def stratum(e): return '>=50' if e >= 50 else ('10-50' if e >= 10 else '<10')
def rd(p):
    d = collections.OrderedDict(); n = None
    for l in open(p):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}
for run in runs:
    d = os.path.join(w, run); cons = rd(os.path.join(d, 'alu.cons.fasta'))
    m = collections.defaultdict(list)
    for l in open(os.path.join(d, 'alu.chunks.tsv')):
        f = l.split(); m[f[1]].append(f[0])
    pur = sum(collections.Counter(truth[x] for x in g).most_common(1)[0][1] for g in m.values()) / S
    st = [l.split('\t') for l in open(os.path.join(d, 'alu.chunkstats.tsv')) if not l.startswith('chunk')]
    agree = sum(float(r[5]) for r in st) / len(st)
    print('== %s: %d rows, %d distinct, chunk purity vs rmsk name %.3f, mean agree %.3f' % (run, len(cons), len(set(cons.values())), pur, agree), flush=True)
    for ref in ('dfam_alu.fa', 'price213.fa'):
        out = os.path.join(d, 'rec_' + ref.split('.')[0])
        if not os.path.exists(out + '.ref_best.tsv'):
            subprocess.run(['python3', os.path.join(here, 'recover.py'), os.path.join(w, ref), os.path.join(d, 'alu.cons.fasta'), out], check=True, stdout=subprocess.DEVNULL)
        rows = [l.rstrip('\n').split('\t') for l in open(out + '.ref_best.tsv')][1:]
        by = collections.defaultdict(list)
        for r in rows: by[stratum(expected(r[0]))].append(float(r[3]))
        print('   %s (%d references):' % (ref, len(rows)))
        for s in ('>=50', '10-50', '<10'):
            v = by.get(s, [])
            if v: print('      expected %-6s n=%3d  recovered at 0.95/0.98/0.99: %3d/%3d/%3d   median best identity %.3f' % (s, len(v), sum(x >= .95 for x in v), sum(x >= .98 for x in v), sum(x >= .99 for x in v), sorted(v)[len(v) // 2]))
        cb = [l.rstrip('\n').split('\t') for l in open(out + '.cons_best.tsv')][1:]
        near = collections.Counter(r[1] for r in cb if float(r[2]) >= 0.95)
        print('      rows with a reference at >= 0.95: %d of %d (nearest references: %d distinct); rows with none: %d' % (sum(near.values()), len(cb), len(near), len(cb) - sum(near.values())), flush=True)
