import subprocess,collections,statistics,sys,glob
rows=[]
for scen in ('young','middle','old'):
  for n in (50,20):
    own_best=[];own_id=[];gap=[];cnt=0
    for seed in (1,2,3,4):
        w='s%d_%s'%(seed,scen)
        cons='%s/sf%d/%s.cons.fasta'%(w,n,scen); fa='%s/%s.fasta'%(w,scen); ch='%s/sf%d/%s.chunks.tsv'%(w,n,scen)
        row={l.split('\t')[0]:l.split('\t')[1] for l in open(ch)}
        r=subprocess.run(['vsearch','--usearch_global',fa,'--db',cons,'--id','0.5','--iddef','1','--strand','both','--maxaccepts','0','--maxrejects','0','--userout','/dev/stdout','--userfields','query+target+id','--quiet','--threads','4'],capture_output=True,text=True)
        ids=collections.defaultdict(dict)
        for l in r.stdout.splitlines():
            q,t,i=l.split('\t'); ids[q][t]=float(i)
        for q,own in row.items():
            d=ids.get(q,{}); 
            if own not in d: own_best.append(0); own_id.append(0.0); continue
            best=max(d.values()); own_best.append(1 if d[own]>=best-1e-9 else 0); own_id.append(d[own]); gap.append(best-d[own])
    own_id_s=sorted(own_id)
    rows.append((scen,n,len(own_id),100*sum(own_best)/len(own_best),statistics.median(own_id),own_id_s[len(own_id_s)//20]))
print('scenario  -n  copies  %copies whose own row is their best row   median id to own row   5th percentile')
for r in rows: print('%-8s %3d %7d %20.1f %28.1f %18.1f'%r)
