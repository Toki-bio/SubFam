#!/usr/bin/env bash
# Ablation on the 3,000 B2 copies: what each part of SubFam contributes. Same copies, N = 50 where chunks exist.
set -e; SUB=$(cd "$(dirname "$0")/../.." && pwd)/SubFam.sh; B=${1:?usage: run_ablation.sh B2_WORKDIR (the make_b2.sh workdir)}; T=8
python3 - "$B/b2_chr19_3000.fa" <<'PY'
import random, sys
seqs=[]; cur=None
for l in open(sys.argv[1]):
    if l.startswith('>'): cur=[l]; seqs.append(cur)
    else: cur.append(l)
random.Random(1).shuffle(seqs); open('shuf.fa','w').write(''.join(''.join(s) for s in seqs))
PY
printf 'variant\trows\tcolumns\tsec\tpurity_rmsk\tpurity_v2021\tmean_agree\tmean_id_min\n' > ablation.tsv
score() { # name dir prefix sec
  python3 - "$1" "$2/$3.chunks.tsv" "$B/b2_truth.tsv" "$B/b2_truth_kramerov.tsv" "$2/$3.chunkstats.tsv" "$2/$3.aln.fasta" "$4" <<'PY'
import sys, collections
name, ch, tr, tk, st, aln, sec = sys.argv[1:8]
m = collections.defaultdict(list)
for l in open(ch): f = l.split(); m[f[1]].append(f[0])
def pur(tf):
    t = dict(l.split()[:2] for l in open(tf)); return sum(collections.Counter(t[x] for x in g).most_common(1)[0][1] for g in m.values()) / sum(len(g) for g in m.values())
rows = [l.split('\t') for l in open(st) if not l.startswith('chunk')]
L = len(''.join(l.strip() for l in open(aln) if not l.startswith('>')).split('>')[0]) if False else None
seq=''; n=0
for l in open(aln):
    if l.startswith('>'): n+=1; 
    elif n==1: seq+=l.strip()
print('%s\t%d\t%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f' % (name, len(m), len(seq), sec, pur(tr), pur(tk), sum(float(r[5]) for r in rows)/len(rows), sum(float(r[7]) for r in rows)/len(rows)))
PY
}
t0=$(date +%s); $SUB -n 50 -t $T -o full -x b2 $B/b2_chr19_3000.fa >/dev/null; score "SubFam (k-mer order, chunks, consensus, alignment)" full b2 $(( $(date +%s)-t0 )) >> ablation.tsv
t0=$(date +%s); $SUB -n 50 -I -t $T -o shuf -x b2 shuf.fa >/dev/null; score "no ordering: shuffled input, same chunks" shuf b2 $(( $(date +%s)-t0 )) >> ablation.tsv
t0=$(date +%s); $SUB -n 50 -I -t $T -o inord -x b2 $B/b2_chr19_3000.fa >/dev/null; score "no ordering: input order (genomic coordinate)" inord b2 $(( $(date +%s)-t0 )) >> ablation.tsv
t0=$(date +%s); $SUB -n 50 -m -t $T -o maf -x b2 $B/b2_chr19_3000.fa >/dev/null; score "MAFFT guide-tree order instead of k-mers" maf b2 $(( $(date +%s)-t0 )) >> ablation.tsv
t0=$(date +%s); $SUB -n 3000 -t $T -o one -x b2 $B/b2_chr19_3000.fa >/dev/null; score "no chunks: one consensus of all copies" one b2 $(( $(date +%s)-t0 )) >> ablation.tsv
# no consensus: align all copies (what the researcher gets without SubFam)
t0=$(date +%s); mafft --thread $T --auto --quiet --reorder $B/b2_chr19_3000.fa > all.aln.fa; sec=$(( $(date +%s)-t0 ))
cols=$(awk '/^>/{if(n++)exit} !/^>/{s=s $0} END{print length(s)}' all.aln.fa)
printf 'no consensus: MAFFT --auto of all 3,000 copies\t3000\t%d\t%d\tNA\tNA\tNA\tNA\n' $cols $sec >> ablation.tsv
cat ablation.tsv

# the last variant (k-mer UPGMA order without the closest-ends flip; ablation only, SUBFAM_NOFLIP=1):
set -e; SUB=$(cd "$(dirname "$0")/../.." && pwd)/SubFam.sh; B=${1:?usage: run_ablation.sh B2_WORKDIR (the make_b2.sh workdir)}; T=8
t0=$(date +%s); SUBFAM_NOFLIP=1 $SUB -n 50 -t $T -o noflip -x b2 $B/b2_chr19_3000.fa >/dev/null; sec=$(( $(date +%s)-t0 ))
python3 - "k-mer UPGMA order without the closest-ends flip" noflip/b2.chunks.tsv $B/b2_truth.tsv $B/b2_truth_kramerov.tsv noflip/b2.chunkstats.tsv noflip/b2.aln.fasta $sec <<'PY' >> ablation.tsv
import sys, collections
name, ch, tr, tk, st, aln, sec = sys.argv[1:8]
m = collections.defaultdict(list)
for l in open(ch): f = l.split(); m[f[1]].append(f[0])
def pur(tf):
    t = dict(l.split()[:2] for l in open(tf)); return sum(collections.Counter(t[x] for x in g).most_common(1)[0][1] for g in m.values()) / sum(len(g) for g in m.values())
rows = [l.split('\t') for l in open(st) if not l.startswith('chunk')]
seq=''; n=0
for l in open(aln):
    if l.startswith('>'): n+=1
    elif n==1: seq+=l.strip()
print('%s\t%d\t%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f' % (name, len(m), len(seq), sec, pur(tr), pur(tk), sum(float(r[5]) for r in rows)/len(rows), sum(float(r[7]) for r in rows)/len(rows)))
PY
