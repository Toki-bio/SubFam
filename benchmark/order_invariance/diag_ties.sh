#!/usr/bin/env bash
# SBATCH -J order_diag
# SBATCH -c 4
# SBATCH --mem=8G
# SBATCH -t 00:40:00
# SBATCH -o /scratch/sk4386/subfam_kmer/diag_%j.out
# Why does the k-mer order depend on input order?  Instrumented copy of kmer_order.c (debug output only; the order it writes is unchanged).
set -euo pipefail
W=/staging/tmp/subfam_order; mkdir -p $W; cd "$W" || exit 1
export TMPDIR=$W
ln -sfn /staging/tmp/subfam_b2_purity/SubFam SubFam
cp -f /staging/tmp/subfam_flanks/f0.fa b2_chr19_3000.fa
for s in $(seq 1 10); do python3 - $s <<'PX'
import sys,random
s=int(sys.argv[1]); r=open('b2_chr19_3000.fa').read().split('>')[1:]; random.Random(100+s).shuffle(r)
open('p%d.fa'%s,'w').write(''.join('>'+x for x in r))
PX
done
git -C SubFam log -1 --format='SubFam commit %h'
python3 - <<'PY'
t=open('SubFam/kmer_order.c').read()
def R(a,b):
    global t
    assert t.count(a)==1,a[:40]; t=t.replace(a,b)
R("static float *CD;","static Rec *RECS;\nstatic uint64_t *HH; static FILE *DBG;\nstatic uint64_t fnv(const char *s){uint64_t h=1469598103934665603ULL;for(;*s;s++){h^=(unsigned char)*s;h*=1099511628211ULL;}return h*0x9E3779B97F4A7C15ULL;}\nstatic float *CD;")
R("    ACTIVE = xmalloc(n);\n","    ACTIVE = xmalloc(n);\n    HH = xmalloc(n * sizeof(uint64_t)); for (size_t i = 0; i < n; i++) HH[i] = fnv(RECS[i].name);\n    DBG = getenv(\"KO_DEBUG\") ? fopen(getenv(\"KO_DEBUG\"), \"w\") : NULL;\n")
R("        long cj = NNJ[ci];\n","""        long cj = NNJ[ci];
        if (DBG) { size_t tg = 0, tr = 0;
            for (size_t i = 0; i < n; i++) if (ACTIVE[i] && NNJ[i] >= 0 && NND[i] == bd) tg++;
            for (size_t j = 0; j < n; j++) if (ACTIVE[j] && (long)j != ci && (double)CD[ci * n + j] == bd) tr++;
            fprintf(DBG, "M %zu %llu %.9g %zu %zu %zu %zu", clen[ci] + clen[cj], (unsigned long long)(HH[ci] + HH[cj]), bd, tg, tr, clen[ci], clen[cj]); fputc(10, DBG); }
""")
R("        int best = 0;\n        for (int t = 1; t < 4; t++) if (opts[t] < opts[best]) best = t;\n","        int best = 0;\n        for (int t = 1; t < 4; t++) if (opts[t] < opts[best]) best = t;\n        if (DBG) { int eq = 0; for (int t = 0; t < 4; t++) if (opts[t] == opts[best]) eq++; fprintf(DBG, \"O %d %zu %zu\", eq, sI, sJ); fputc(10, DBG); }\n")
R("        ACTIVE[cj] = 0;\n","        ACTIVE[cj] = 0; HH[ci] += HH[cj];\n")
R("Rec *recs = read_fasta(argv[1], &n);","Rec *recs = read_fasta(argv[1], &n); RECS = recs;")
open('kmer_order_dbg.c','w').write(t)
PY
gcc -O2 -o kmer_order_dbg kmer_order_dbg.c -lm
gcc -O2 -o kmer_order_ref SubFam/kmer_order.c -lm
for f in b2_chr19_3000 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10; do
  KO_DEBUG=dbg_$f.tsv ./kmer_order_dbg $f.fa ord_$f.fa 6 0 0
done
# the debug build must write the same order as the release code
./kmer_order_ref b2_chr19_3000.fa ord_ref.fa 6 0 0; cmp ord_ref.fa ord_b2_chr19_3000.fa && echo "debug build order == release order"
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
C0,gt,rt,ot,zero=clades('dbg_%s.tsv'%base)
print('merges %d: global min tied %d, tied inside the chosen row %d, zero-distance merges %d, flip options tied %d'%(sum(C0.values()),gt,rt,zero,ot))
o0=names('ord_%s.fa'%base)
def adj(o): return set(frozenset(p) for p in zip(o,o[1:]))
A0=adj(o0)
print('run\tclades shared with clean run\tleaf order identical up to reversal\tneighbour pairs shared\ttied global\ttied row\ttied flip')
for f in ['p%d'%i for i in range(1,11)]:
    C,gt,rt,ot,zero=clades('dbg_%s.tsv'%f); o=names('ord_%s.fa'%f)
    sh=sum((C&C0).values())/sum(C0.values()); same=(o==o0 or o==o0[::-1])
    print('%s\t%.4f\t%s\t%.3f\t%d\t%d\t%d'%(f,sh,same,len(A0&adj(o))/len(A0),gt,rt,ot))
PY
echo ALL_DONE
