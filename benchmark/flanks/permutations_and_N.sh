#!/usr/bin/env bash
#SBATCH -J subfam_kmer
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 02:00:00
#SBATCH -o /scratch/sk4386/subfam_kmer/slurm_%j.out
# Monsoon (Slurm, 4 threads), 2026-11-01, SubFam 0669447. (1) does breaking k-mers at N change anything?  (2) 10 input permutations of the 3,000 B2 copies. Needs cmp_funcs.py = the functions of compare.py in flanks.sh.
set -euo pipefail
W=/scratch/sk4386/subfam_kmer; mkdir -p $W/dl; cd "$W" || exit 1
export TMPDIR=$W/tmp; mkdir -p $TMPDIR
source ~/.bashrc 2>/dev/null || true
export PATH=/home/sk4386/.conda/envs/bioinfo/bin:$PATH
python3 -c "import Bio" 2>/dev/null || python3 -m pip install -q biopython
[ -d SubFam ] || git clone -q https://github.com/Toki-bio/SubFam.git
TOP=$W/SubFam; git -C $TOP log -1 --format='SubFam commit %h'
[ -s dl/rmsk.txt.gz ] || curl -sSf -o dl/rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz
[ -s dl/chr19.fa.gz ] || curl -sSf -o dl/chr19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr19.fa.gz
python3 -I $TOP/benchmark/line_real/extract_l1.py dl/rmsk.txt.gz dl/chr19.fa.gz b2_chr19.fa 170 '^(B2_Mm1a|B2_Mm1t|B2_Mm2|B3|B3A)$'
python3 -I $TOP/benchmark/sine_b2/subsample.py b2_chr19.fa 3000 1 b2_chr19_3000.fa b2_truth.tsv
grep -c '>' b2_chr19_3000.fa
# fixed variant: the k-mer window restarts after a letter that is not A/C/G/T/U
rm -rf SubFam_fix; cp -r $TOP SubFam_fix
python3 - <<'PY'
p='SubFam_fix/kmer_order.c'; t=open(p).read()
a="        if (c < 0) continue;\n        h = ((h << 2)"; assert t.count(a)==1
open(p,'w').write(t.replace(a,"        if (c < 0) { have = 0; continue; }\n        h = ((h << 2)"))
PY
# N-bearing input: 1,000 random copies get one run of 1-10 N in the interior (same Ns for both variants)
python3 - <<'PY'
import random
r=random.Random(7); recs=open('b2_chr19_3000.fa').read().split('>')[1:]
out=[]; idx=set(r.sample(range(len(recs)),1000)); n=0
for i,x in enumerate(recs):
    h,*s=x.split('\n'); s=''.join(s)
    if i in idx:
        L=r.randint(1,10); a=r.randint(20,len(s)-20-L); s=s[:a]+'N'*L+s[a+L:]; n+=1
    out.append('>%s\n%s\n'%(h,s))
open('b2_withN.fa','w').write(''.join(out)); print('copies with an N run',n)
PY
run(){ # tree out prefix fasta
  ( cd $W; $1/SubFam.sh -n 50 -t 4 -o $2 -x $2 $3 > $2.log 2>&1 ); }
run $TOP        o_clean b2_chr19_3000.fa
run $TOP        o_N     b2_withN.fa
run $W/SubFam_fix f_N   b2_withN.fa
run $W/SubFam_fix f_clean b2_chr19_3000.fa
# permutations of the clean input (old code)
for s in $(seq 1 10); do
  python3 - $s <<'PY'
import sys,random
s=int(sys.argv[1]); recs=open('b2_chr19_3000.fa').read().split('>')[1:]; random.Random(100+s).shuffle(recs)
open('p%d.fa'%s,'w').write(''.join('>'+x for x in recs))
PY
  run $TOP p$s p$s.fa
done
cp $TOP/benchmark/sine_b2/kramerov2021_B2_consensus.fa . 2>/dev/null || true
python3 - <<'PY'
import statistics as st
exec(open('cmp_funcs.py').read())
P=lambda r:part('%s/%s.chunks.tsv'%(r,r)); C=lambda r:fa('%s/%s.cons.fasta'%(r,r))
print('--- k-mers across N: chunk partitions (pairs of the first run kept in the same chunk by the second)')
for a,b in (('o_clean','o_N'),('o_clean','f_N'),('o_N','f_N'),('o_clean','f_clean')):
    print(a,'->',b,'pairs kept %.3f  ARI %.3f'%(kept_pairs(P(a),P(b)),ari(P(a),P(b))))
r0=C('o_clean')
def rec(run):
    r=C(run); b=sorted(best(s,r) for s in r0.values()); return sum(x>=0.95 for x in b), b[30]
print('--- original rows (clean old code) found in:')
for run in ('o_N','f_N','f_clean'): print(run,'%d/60 at >=0.95, median best id %.3f'%rec(run))
print('--- 10 input permutations of the clean input (old code)')
res=[]; kp=[]
for s in range(1,11):
    f,m=rec('p%d'%s); res.append(f); kp.append(kept_pairs(P('o_clean'),P('p%d'%s)))
    print('p%d recurring %d/60 median %.3f pairs kept %.3f'%(s,f,m,kp[-1]))
print('recurring rows: min %d median %d max %d | pairs kept: min %.3f median %.3f max %.3f'%(min(res),st.median(res),max(res),min(kp),st.median(kp),max(kp)))
PY
echo ALL_DONE
