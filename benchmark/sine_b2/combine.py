#!/usr/bin/env python3
"""combine.py OUT.aln.fasta COPIES.fa CHUNKS.tsv CONS.fa TRUTH1.tsv NAME1 TRUTH2.tsv NAME2 REFS.fa... [--threads T --tmp DIR]

One alignment for manual inspection: every published consensus (REFS, in the given order) at the
top, then the SubFam chunk consensuses in chunk order (= k-mer guide-tree order), each named
  <chunk>|<NAME1>:<majority label>_p<purity>|<NAME2>:<majority label>_p<purity>
so that two classifications of the same copies can be read side by side. MAFFT --auto."""
import collections, os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'examples'))
from build_examples import rd, mafft

args = sys.argv[1:]; threads, tmp = 4, None
if '--threads' in args: i = args.index('--threads'); threads = int(args[i + 1]); del args[i:i + 2]
if '--tmp' in args: i = args.index('--tmp'); tmp = args[i + 1]; del args[i:i + 2]
out, copies, chunks, cons, t1, n1, t2, n2 = args[:8]; refs = args[8:]
tmp = tmp or tempfile.mkdtemp(); os.makedirs(tmp, exist_ok=True)
truth = [dict(l.split()[:2] for l in open(t) if l.strip()) for t in (t1, t2)]
members = collections.defaultdict(list)
for l in open(chunks):
    f = l.split()
    if len(f) >= 2: members[f[1]].append(f[0])
rows = collections.OrderedDict()
for r in refs:
    for k, v in rd(r).items(): rows['PUBLISHED_' + k] = v.upper()
for ch, seq in rd(cons).items():
    tag = []
    for name, tr in ((n1, truth[0]), (n2, truth[1])):
        labs = [tr[m] for m in members[ch] if m in tr]
        if labs:
            sf, n = collections.Counter(labs).most_common(1)[0]; tag.append('%s:%s_p%.2f' % (name, sf, n / len(labs)))
        else: tag.append(name + ':none')
    rows['%s|%s' % (ch, '|'.join(tag))] = seq
aln = mafft(rows, threads, tmp, 'combined', ('--auto',))
with open(out, 'w') as o:
    for k in rows: o.write('>%s\n%s\n' % (k, aln[k]))
print('wrote', out, len(rows), 'rows', file=sys.stderr)
