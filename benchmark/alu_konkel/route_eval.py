import collections, json, subprocess, os, sys
EV='/home/user/SubFam/benchmark/evaluate.py'
SUB='/home/user/SubFam/SubFam.sh'
lab={}
for l in open('../alutest/loci_labels.tsv'):
    f=l.rstrip('\n').split('\t')
    if f[0]!='acc': lab[f[0]]=f[6]
def rd(p):
    n=[];s=[]
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n.append(l[1:].split()[0]);s.append([])
        else: s[-1].append(l.strip())
    return dict(zip(n,[''.join(x) for x in s]))
bodies=rd('bodies295.fa')                       # id like KT305466|AluY
def acc(i): return i.split('|')[0]
ids=[i for i in bodies if lab[acc(i)]!='ambig']
truth={i:lab[acc(i)] for i in ids}
os.makedirs('ev',exist_ok=True)
open('ev/truth.tsv','w').write(''.join('%s\t%s\n'%(i,t) for i,t in truth.items()))
# masters: the three lineages present among the labels
P=rd('../price213.fa'); m=open('ev/masters.fa','w')
for k,v in P.items():
    if k.split('|')[0] in ('AluY','AluYa5','AluYb8'): m.write('>%s\n%s\n'%(k.split('|')[0],v))
m.close()
def consensus(copy_ids,tag):
    f='ev/%s.fa'%tag; open(f,'w').write(''.join('>%s\n%s\n'%(i,bodies[i]) for i in copy_ids))
    if len(copy_ids)==1: return bodies[copy_ids[0]]
    subprocess.run([SUB,'-n',str(len(copy_ids)+1),'-t','4','-o','ev/'+tag+'_o','-x','c',f],capture_output=True)
    r=rd('ev/%s_o/c.cons.fasta'%tag); return list(r.values())[0]
def score(name,groups):                         # groups: gid -> list of copy ids (non-ambig only)
    groups={g:v for g,v in groups.items() if v}
    fa=open('ev/reps_%s.fa'%name.replace(' ','_'),'w'); mem=open('ev/mem_%s.tsv'%name.replace(' ','_'),'w')
    for g,v in groups.items():
        c=CONS[name][g] if name in CONS and g in CONS[name] else consensus(v,name.replace(' ','_')+'_'+g)
        fa.write('>%s\n%s\n'%(g,c))
        for i in v: mem.write('%s\t%s\n'%(i,g))
    fa.close(); mem.close()
    r=subprocess.run(['python3','-I',EV,name,fa.name,mem.name,'ev/truth.tsv','ev/masters.fa'],capture_output=True,text=True)
    print(r.stdout.strip() or r.stderr[-300:])
CONS={}
def chunks(path):
    g=collections.defaultdict(list)
    for l in open(path):
        i,c,s=l.rstrip('\n').split('\t')
        if i in truth: g[c].append(i)
    return g
print('method\treps\treps_ge10\tsingleton_frac\tpurity\tmedian_id\tsubfam_recovered(of 3)')
score('SubFam_n20',chunks('sf20/sf.chunks.tsv'))
score('SubFam_n5',chunks('sf5/sf.chunks.tsv'))
# peel route on n=5 chunk consensuses
ch=chunks('sf5/sf.chunks.tsv'); pj=json.load(open('peel/peel_features.json'))
grp={}; used=set()
for k,p in enumerate(pj['peeled']):
    grp['peel%d'%(k+1)]=[i for c in p['members'] for i in ch.get(c,[])]; used|=set(p['members'])
grp['residue']=[i for c in pj['residue'] for i in ch.get(c,[])]
score('SubFam_n5_peel',grp)
print('peel groups (copies, labels):',{g:dict(collections.Counter(truth[i] for i in v)) for g,v in grp.items()})
# COSEG
rows=[l.rstrip('\n') for l in open('../cosegrun/m50/konkel.seqs')]
names=[l.strip() for l in open('../cosegrun/konkel.names')]
for m in (50,5):
    asg=[l.split()[-1] for l in open('../cosegrun/m%d/konkel.seqs.assign'%m)]
    g=collections.defaultdict(list)
    for n,a in zip(names,asg):
        if n in truth: g['g'+a].append(n)
    cs=rd('../cosegrun/cons_m%d.fa'%m)
    CONS['COSEG_m%d'%m]={'g'+k.split('_g')[1].split('_')[0]:v for k,v in cs.items()}
    score('COSEG_m%d'%m,g)
# VSEARCH
for idv in ('0.90','0.95','0.98','0.99'):
    subprocess.run(['vsearch','--cluster_fast','bodies295.fa','--id',idv,'--uc','ev/vs.uc','--centroids','ev/vs_c.fa','--quiet','--threads','4'])
    cent={}; g=collections.defaultdict(list)
    for l in open('ev/vs.uc'):
        f=l.split('\t')
        if f[0]=='S': cent[f[8]]=f[8]; 
    for l in open('ev/vs.uc'):
        f=l.rstrip('\n').split('\t')
        if f[0] in 'SH':
            c=f[9] if f[0]=='H' else f[8]
            if f[8] in truth: g[c].append(f[8])
    CONS['VSEARCH_%s'%idv]={c:bodies[c] for c in g}
    score('VSEARCH_%s'%idv,g)
