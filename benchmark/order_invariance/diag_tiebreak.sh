#!/usr/bin/env bash
# Same diagnostic with a content-based tie-break (experiment): does the order stop depending on the input order?
set -euo pipefail
cd /staging/tmp/subfam_order || exit 1
cp kmer_order_dbg.c kmer_order_base.c
python3 tb_patch.py kmer_order_base.c kmer_order_tb.c
gcc -O2 -o kmer_order_tb kmer_order_tb.c -lm
for f in b2_chr19_3000 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10; do
  KO_DEBUG=tb_dbg_$f.tsv ./kmer_order_tb $f.fa tb_ord_$f.fa 6 0 0
done
python3 - <<'PY'
import collections
def names(f): return [l[1:].split()[0] for l in open(f) if l.startswith('>')]
def seqs(f):
    d={};n=None
    for l in open(f):
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l.strip().upper())
    return {k:''.join(v) for k,v in d.items()}
base='b2_chr19_3000'; S=seqs(base+'.fa')
dup=collections.Counter(S.values()); print('identical-sequence copies: %d copies share a sequence with another (%d distinct sequences of %d)'%(sum(v for v in dup.values() if v>1),len(dup),len(S)))
def clades(f):
    c=collections.Counter(); ties=0; gt=0; rt=0; ot=0; zero=0
    for l in open(f):
        x=l.split()
        if x[0]=='M':
            c[(int(x[1]),int(x[2]))]+=1; g,r=int(x[4]),int(x[5])
            if g>1: gt+=1
            if r>1: rt+=1
            if float(x[3])==0: zero+=1
        else:
            if int(x[1])>1: ot+=1
    return c,gt,rt,ot,zero
C0,gt,rt,ot,zero=clades('tb_dbg_%s.tsv'%base)
print('merges %d: global min tied %d, tied inside the chosen row %d, zero-distance merges %d, flip options tied %d'%(sum(C0.values()),gt,rt,zero,ot))
o0=names('tb_ord_%s.fa'%base)
def adj(o): return set(frozenset(p) for p in zip(o,o[1:]))
A0=adj(o0)
print('run\tclades shared with clean run\tleaf order identical up to reversal\tneighbour pairs shared\ttied global\ttied row\ttied flip')
for f in ['p%d'%i for i in range(1,11)]:
    C,gt,rt,ot,zero=clades('tb_dbg_%s.tsv'%f); o=names('tb_ord_%s.fa'%f)
    sh=sum((C&C0).values())/sum(C0.values()); same=(o==o0 or o==o0[::-1])
    print('%s\t%.4f\t%s\t%.3f\t%d\t%d\t%d'%(f,sh,same,len(A0&adj(o))/len(A0),gt,rt,ot))
PY
echo ALL_DONE
