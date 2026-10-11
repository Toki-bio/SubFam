#!/usr/bin/env python3
"""extract_alu.py rmsk.txt.gz chrN.fa.gz OUT.fa [MINLEN=250] [MAXLEN=330]
Full-length Alu copies of one hg38 chromosome from the UCSC rmsk table (repFamily == 'Alu', MINLEN <= length <= MAXLEN, no N),
on the element strand, named <repName>_<chr>_<start>_<end>. Python-only counterpart of extract_alu.sh (no bedtools/seqkit).
The rmsk repName is kept in the id as an external label; SubFam does not use it."""
import sys, gzip
rm, fa, out = sys.argv[1:4]; mn = int(sys.argv[4]) if len(sys.argv) > 4 else 250; mx = int(sys.argv[5]) if len(sys.argv) > 5 else 330
name = None; buf = []
with gzip.open(fa, 'rt') as f:
    for l in f:
        if l.startswith('>'):
            if name: break
            name = l[1:].split()[0]
        else: buf.append(l.strip())
chrom = name; seq = ''.join(buf).upper(); del buf
comp = str.maketrans('ACGTN', 'TGCAN'); n = 0; byname = {}
with gzip.open(rm, 'rt') as f, open(out, 'w') as o:
    for l in f:
        c = l.rstrip('\n').split('\t')
        if c[5] != chrom or c[12] != 'Alu': continue
        s, e = int(c[6]), int(c[7])
        if not mn <= e - s <= mx: continue
        t = seq[s:e]
        if 'N' in t: continue
        if c[9] == '-': t = t.translate(comp)[::-1]
        o.write('>%s_%s_%d_%d\n%s\n' % (c[10], chrom, s, e, t)); n += 1; byname[c[10]] = byname.get(c[10], 0) + 1
print(chrom, n, 'copies;', len(byname), 'repNames; top:', sorted(byname.items(), key=lambda x: -x[1])[:12], file=sys.stderr)
