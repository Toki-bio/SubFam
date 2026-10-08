import glob,re,sys,collections
from Bio import Align
from Bio.Seq import Seq
D,PRICE,OUT=sys.argv[1:4]
cons={}
for l in open(PRICE):
    m=re.match(r'subfamily\s+(\d+)\s+(\S+)\s+(\d+)\s+\S+(?:\s+\([^)]*\))?\s+([acgtn]+)\s*$',l)
    if m: cons[m.group(2)]=m.group(4).upper()
cand=['AluY','AluYa5','AluYb8','AluYb9']
al=Align.PairwiseAligner(); al.mode='local'; al.match_score=2; al.mismatch_score=-3; al.open_gap_score=-5; al.extend_gap_score=-2
recs=[]
for f in sorted(glob.glob(D+'/KT*.gb')):
    t=open(f).read()
    acc=re.search(r'^ACCESSION\s+(\S+)',t,re.M).group(1)
    seq=''.join(re.findall(r'^\s+\d+\s+([a-z ]+)$',t.split('ORIGIN')[1],re.M)).replace(' ','').upper()
    d=' '.join(re.search(r'^DEFINITION\s+(.*?)\n(?=\S)',t,re.M|re.S).group(1).split())
    recs.append((acc,d,seq))
rows=[];fa=open(OUT+'/alu_bodies.fa','w')
for acc,d,seq in recs:
    best=None
    for strand,s in (('+',seq),('-',str(Seq(seq).reverse_complement()))):
        a=al.align(cons['AluY'],s)[0]   # query=consensus, target=locus
        sc=a.score
        if best is None or sc>best[0]: best=(sc,strand,s,a)
    sc,strand,s,a=best
    ts=a.aligned[1][0][0]; te=a.aligned[1][-1][1]; qs=a.aligned[0][0][0]; qe=a.aligned[0][-1][1]
    body=s[ts:te]
    cov=(qe-qs)/len(cons['AluY'])
    sc2={}
    for c in cand:
        sc2[c]=al.align(cons[c],body)[0].score
    order=sorted(sc2,key=lambda c:-sc2[c])
    lab=order[0] if sc2[order[0]]>sc2[order[1]] else 'ambig'
    rows.append((acc,strand,len(seq),ts,te,round(cov,2),lab,sc2['AluY'],sc2['AluYa5'],sc2['AluYb8'],sc2['AluYb9'],d))
    if cov>=0.9 and 250<=len(body)<=330:
        fa.write('>%s|%s\n%s\n'%(acc,lab,body))
with open(OUT+'/loci_labels.tsv','w') as o:
    o.write('acc\tstrand\tlocus_len\tbody_start\tbody_end\tcov_vs_AluY\tlabel\tsAluY\tsYa5\tsYb8\tsYb9\tdef\n')
    for r in rows: o.write('\t'.join(map(str,r))+'\n')
print(len(rows),'loci;',collections.Counter(r[6] for r in rows)); print('full-length kept:',sum(1 for r in rows if r[5]>=0.9),'strand',collections.Counter(r[1] for r in rows))
