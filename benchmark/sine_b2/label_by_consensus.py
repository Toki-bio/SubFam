#!/usr/bin/env python3
"""label_by_consensus.py COPIES.fa CONSENSUS.fa OUT.tsv [RMSK_TRUTH.tsv] [--threads T]

Best-hit label of every copy against a set of consensus sequences, the way Borodulina et al. 2025
assigned B2 copies to subfamilies (best score to a subfamily consensus): vsearch --usearch_global,
identity over the alignment (--iddef 2), threshold 0.5, all hits kept and the best one(s) taken
(a tie is 'tie_<a>_<b>'). The word prefilter is loosened (--wordlength 5 --minwordmatches 3) so
that copies 30 % diverged from every consensus still get a hit. Writes id<TAB>label; with a second
truth file prints the cross-table of the two labellings and the median identity per label."""
import collections, subprocess, sys

args = sys.argv[1:]; threads = 4
if '--threads' in args: i = args.index('--threads'); threads = int(args[i + 1]); del args[i:i + 2]
copies, cons, out = args[:3]; rmsk = args[3] if len(args) > 3 else None
hits_file = out + '.hits'
subprocess.run(['vsearch', '--usearch_global', copies, '--db', cons, '--id', '0.5', '--iddef', '2', '--strand', 'plus',
                '--maxaccepts', '0', '--maxrejects', '0', '--top_hits_only', '--wordlength', '5', '--minwordmatches', '3',
                '--userout', hits_file, '--userfields', 'query+target+id', '--quiet', '--threads', str(threads)], check=True)
names = [l[1:].split()[0] for l in open(copies) if l.startswith('>')]
hits = collections.defaultdict(list)
for l in open(hits_file):
    q, t, i = l.split(); hits[q].append((t.split('_', 1)[1] if '_' in t else t, float(i)))
lab, ids = {}, {}
for q, h in hits.items():
    best = max(x[1] for x in h); tops = sorted({t for t, i in h if i == best})
    lab[q] = tops[0] if len(tops) == 1 else 'tie_' + '_'.join(tops); ids[q] = best
with open(out, 'w') as o:
    for q in names: o.write('%s\t%s\n' % (q, lab.get(q, 'none')))
print('copies with a hit:', len(lab), 'of', len(names), '; ties:', sum(v.startswith('tie') for v in lab.values()), file=sys.stderr)
byk = collections.defaultdict(list)
for q in ids: byk[lab[q]].append(ids[q])
for k in sorted(byk): v = sorted(byk[k]); print('%s\t%d copies\tmedian identity to consensus %.1f' % (k, len(v), v[len(v) // 2]), file=sys.stderr)
if rmsk:
    truth = dict(l.split()[:2] for l in open(rmsk) if l.strip())
    ct = collections.Counter((truth[q], lab.get(q, 'none')) for q in names if q in truth)
    rm = sorted({k[0] for k in ct}); kr = sorted({k[1] for k in ct})
    print('rmsk \\ consensus', *kr, sep='\t')
    for r in rm: print(r, *[ct.get((r, k), 0) for k in kr], sep='\t')
