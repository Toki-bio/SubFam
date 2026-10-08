#!/usr/bin/env python3
"""Cut RepeatMasker L1 copies of given subfamilies out of an hg38 chromosome (UCSC rmsk table + chromosome fasta).
usage: extract_l1.py rmsk.txt.gz chrN.fa.gz OUT.fa [MINLEN] [REGEX]     defaults: 1000, '^L1PA[2-8]$'
Copies are on the + strand of the element (reverse-complemented where rmsk says '-'), upper case,
named <repName>_<chr>_<start>_<end> (0-based, as in rmsk). The rmsk call is kept in the name: it is an external label, not used by SubFam."""
import sys, gzip, re
rm, fa, out = sys.argv[1:4]; minlen = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
rx = re.compile(sys.argv[5] if len(sys.argv) > 5 else r'^L1PA[2-8]$')
name = None; buf = []
with gzip.open(fa, 'rt') as f:
    for l in f:
        if l.startswith('>'):
            if name: break
            name = l[1:].split()[0]
        else: buf.append(l.strip())
chrom = name; seq = ''.join(buf).upper(); del buf
comp = str.maketrans('ACGTN', 'TGCAN')
n = 0
with gzip.open(rm, 'rt') as f, open(out, 'w') as o:
    for l in f:
        c = l.rstrip('\n').split('\t')
        if c[5] != chrom or not rx.match(c[10]): continue
        s, e = int(c[6]), int(c[7])
        if e - s < minlen: continue
        t = seq[s:e]
        if 'N' in t: continue
        if c[9] == '-': t = t.translate(comp)[::-1]
        o.write('>%s_%s_%d_%d\n%s\n' % (c[10], chrom, s, e, t)); n += 1
print(chrom, n, 'copies', file=sys.stderr)
