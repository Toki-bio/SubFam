# Converts the Alu bodies (prep.py output) to COSEG input: aligns each body to a consensus, writes NAME.seqs / NAME.ins / NAME.names.
# usage: to_coseg.py ALU.cons alu_bodies.fa OUTPREFIX   (ALU.cons from github.com/rmhubley/coseg; copies missing >5 bases at an end are dropped, as in Price 2004)
import sys,re
from Bio import Align
cons_fa,bodies,out=sys.argv[1:4]
cons=''.join(l.strip() for l in open(cons_fa) if not l.startswith('>'))
L=len(cons)
al=Align.PairwiseAligner(); al.mode='global'; al.match_score=2; al.mismatch_score=-3
al.open_gap_score=-5; al.extend_gap_score=-2; al.target_end_gap_score=0; al.query_end_gap_score=0
names=[];seqs=[];ins=[];kept=0;drop=0
recs=[]
cur=None
for l in open(bodies):
    if l.startswith('>'): recs.append([l[1:].strip(),'']); 
    else: recs[-1][1]+=l.strip()
fs=open(out+'.seqs','w'); fi=open(out+'.ins','w'); fn=open(out+'.names','w')
for name,s in recs:
    a=al.align(cons.upper(),s.upper())[0]  # target=cons, query=body
    row=['-']*L; insd={}
    (tb,qb)=a.aligned
    # walk blocks
    tpos=0
    prev_t=None;prev_q=None
    for (t0,t1),(q0,q1) in zip(tb,qb):
        for k in range(t1-t0): row[t0+k]=s[q0+k].lower()
        if prev_t is not None and q0>prev_q:       # insertion in body between blocks
            insd[prev_t]=s[prev_q:q0].lower()
        prev_t,prev_q=t1,q1
    # trim rule: copies missing >5 bases at either end are dropped
    lead=len(row)-len(''.join(row).lstrip('-')); trail=len(row)-len(''.join(row).rstrip('-'))
    if lead>5 or trail>5: drop+=1; continue
    kept+=1
    line=[];insl=[]
    for i,c in enumerate(row):
        line.append(c)
        if (i+1) in insd: line.append('+'); insl.append('%d:%s'%(i+1,insd[i+1]))
    fs.write(''.join(line)+'\n'); fi.write(' '.join(insl)+'\n'); fn.write(name+'\n')
print('kept',kept,'dropped (>5 bases missing at an end)',drop)
