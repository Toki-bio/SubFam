#!/usr/bin/env python3
"""Build a viewer alignment: SubFam consensus rows grouped by the subfamily they mostly come from,
each group headed by the consensus of its own member copies (the best consensus a method could reach,
CpG decay and all), not the master. Reference rows are named TRUE_<sf>; SubFam rows <chunk>|<sf>|p<purity>.

usage: build_examples.py OUT.aln.fasta COPIES.fa TRUTH.tsv CONS.fa CHUNKS.tsv [--minspan N] [--threads T] [--tmp DIR]
  TRUTH.tsv   copy<TAB>subfamily       CHUNKS.tsv  copy<TAB>chunk<TAB>strand (SubFam's *.chunks.tsv)
"""
import random, sys, subprocess, collections, os, argparse, tempfile

def rd(f):
    d = collections.OrderedDict(); n = None
    for l in open(f):
        l = l.rstrip()
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return collections.OrderedDict((k, ''.join(v)) for k, v in d.items())

def mafft(seqs, threads, tmp, tag, extra=()):
    p = os.path.join(tmp, tag + '.in.fa')
    with open(p, 'w') as o:
        for k, v in seqs.items(): o.write('>%s\n%s\n' % (k, v))
    r = subprocess.run(['mafft', '--quiet', '--thread', str(threads), *extra, p], capture_output=True, text=True, check=True)
    out = collections.OrderedDict(); n = None
    for l in r.stdout.splitlines():
        if l.startswith('>'): n = l[1:].split()[0]; out[n] = []
        elif n: out[n].append(l.strip())
    return collections.OrderedDict((k, ''.join(v).upper()) for k, v in out.items())

def plurality(aln, frac=0.36, mincov=3):
    """coverage-relative plurality (SubFam -c): candidates are the sequences spanning the column"""
    rows = list(aln.values()); L = len(rows[0]); spans = []
    for s in rows:
        a = len(s) - len(s.lstrip('-')); b = len(s.rstrip('-')); spans.append((a, b))
    out = []
    for i in range(L):
        cnt = collections.Counter(); ncov = 0
        for s, (a, b) in zip(rows, spans):
            if a <= i < b:
                ncov += 1; cnt[s[i]] += 1
        if ncov < min(mincov, len(rows)): continue
        need = int(frac * ncov + 1e-9); need += (frac * ncov - need > 1e-9)   # ceil, safe for 0.07 * 100
        ch, c = max(((k, v) for k, v in cnt.items()), key=lambda kv: (kv[1], kv[0] != '-'))
        if ch != '-' and c >= need: out.append(ch)
    return ''.join(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('copies'); ap.add_argument('truth'); ap.add_argument('cons'); ap.add_argument('chunks')
    ap.add_argument('--minspan', type=int, default=0, help='use only copies at least this long for the reference consensus')
    ap.add_argument('--sample', type=int, default=0, help='build each member consensus from at most this many randomly chosen copies (seed 1); 0 = all')
    ap.add_argument('--threads', type=int, default=4); ap.add_argument('--tmp', default=None)
    ap.add_argument('--prefix', default='TRUE')
    ap.add_argument('--refs', default=None, help='FASTA of published consensuses named <anything>_<subfamily>; used instead of member consensuses')
    ap.add_argument('--both', action='store_true', help='with --refs: build the member consensuses too, and put each published consensus right after the head of its group')
    ap.add_argument('--extra-heads', default=None, help='FASTA of further published consensuses, added as REFERENCE_<name> rows at the top, in file order; not tied to any group, ignored by check_rows.py')
    a = ap.parse_args()
    tmp = a.tmp or tempfile.mkdtemp(); os.makedirs(tmp, exist_ok=True)
    copies = rd(a.copies); truth = dict(l.split()[:2] for l in open(a.truth) if l.strip())
    cons = rd(a.cons); chunk_of = {}
    members = collections.defaultdict(list)
    for l in open(a.chunks):
        f = l.split()
        if len(f) >= 2 and f[0] in truth: chunk_of[f[0]] = f[1]; members[f[1]].append(truth[f[0]])
    rowinfo = {}
    for ch, labs in members.items():
        c = collections.Counter(labs); sf, n = c.most_common(1)[0]; rowinfo[ch] = (sf, n / len(labs))
    sfs = sorted({v[0] for v in rowinfo.values()})
    refs = collections.OrderedDict()
    pubname, pub = {}, {}
    for k, v in (rd(a.refs).items() if a.refs else []):
        sf = next((x for x in sfs if k == x or k.endswith('_' + x)), None)
        if sf: pubname[sf] = k; pub[sf] = v.upper()
        else: pubname[k] = k; pub[k] = v.upper()
    for sf in sfs:
        if sf in pub and not a.both: refs[sf] = pub[sf]; continue
        if a.refs and not a.both: continue
        mem = {k: copies[k] for k in copies if truth.get(k) == sf and len(copies[k]) >= a.minspan}
        if len(mem) < 3: continue
        if a.sample and len(mem) > a.sample:
            keep = sorted(random.Random(1).sample(sorted(mem), a.sample)); mem = {k: mem[k] for k in keep}
        aln = mafft(mem, a.threads, tmp, 'ref_' + sf, ('--auto',))
        refs[sf] = plurality(aln)
        print(sf, len(mem), 'copies ->', len(refs[sf]), 'bp', file=sys.stderr)
    rows = collections.OrderedDict()
    for k, v in (rd(a.extra_heads).items() if a.extra_heads else []): rows['REFERENCE_' + k] = v.upper()
    for sf in sfs:
        if sf in refs: rows[('PUBLISHED_%s' % pubname[sf]) if sf in pub and not a.both else ('%s_%s_consensus_of_%s_copies' % (a.prefix, sf, sf))] = refs[sf]
        if a.both and sf in pub: rows['PUBLISHED_%s' % pubname[sf]] = pub[sf]
        for ch in cons:
            if ch in rowinfo and rowinfo[ch][0] == sf:
                rows['%s|%s|p%.2f' % (ch, sf, rowinfo[ch][1])] = cons[ch]
    for sf in pub:
        if sf not in sfs: rows['PUBLISHED_%s_no_group_has_this_majority' % pubname[sf]] = pub[sf]
    for ch in cons:
        if ch not in rowinfo: rows[ch + '|none'] = cons[ch]
    aln = mafft(rows, a.threads, tmp, 'final', ('--auto',))
    with open(a.out, 'w') as o:
        for k in rows: o.write('>%s\n%s\n' % (k, aln[k]))
    print('wrote', a.out, len(rows), 'rows', file=sys.stderr)

if __name__ == '__main__':
    main()
