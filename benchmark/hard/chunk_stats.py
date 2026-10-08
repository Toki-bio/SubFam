#!/usr/bin/env python3
"""Per-chunk quality statistics from the kept chunk alignments of a SubFam run (-K).

usage: chunk_stats.py WORKDIR PREFIX [TRUTH.tsv]
WORKDIR = the .subfam_* directory kept by -K (chunk alignments PREFIX_NN.fasta.aln); output goes to stdout as a table.
Per chunk:
  n         sequences
  id_med    median identity of the members to the chunk's plurality sequence over core columns (columns where >= 50 % of members have a base)
  agree     mean fraction of covering members that carry the plurality base in a core column (how tight the alignment is)
  ragged    fraction of core columns where >= 20 % of the members have a gap (alignment quality)
  f2, best  2-means split of the members on variable columns (minority >= 15 %): size of the smaller part, and the strongest single column that
            separates the two parts (|difference of minority-base frequency|, as in SINE-discriminator/pairsurvey.py)
  outl      fraction of members whose identity is below median - max(3 MAD, 3 points)
  type      T1 one lineage (no split with f2 >= 0.2 and best >= 0.7, outl < 0.10), T2 two lineages (split), TJ one lineage plus junk (outl >= 0.10)
With TRUTH.tsv (copy id -> subfamily, simulations only) also prints the true composition class: pure (top share >= 0.9), two (top two each >= 0.2 and together >= 0.9), other."""
import collections, glob, os, statistics, sys
import numpy as np

def read_aln(p):
    names, seqs = [], []
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'): names.append(l[1:].split()[0]); seqs.append([])
        elif seqs: seqs[-1].append(l.strip().upper())
    return names, np.array([list(''.join(s)) for s in seqs])

def stats(A):
    n, L = A.shape
    nongap = A != '-'
    core = nongap.mean(axis=0) >= 0.5
    if core.sum() < 20: return None
    B = A[:, core]; G = nongap[:, core]
    cons = []; agree = []; ragged = 0
    for j in range(B.shape[1]):
        c = collections.Counter(B[G[:, j], j]); b, k = c.most_common(1)[0]
        cons.append(b); agree.append(k / G[:, j].sum())
        if (~G[:, j]).mean() >= 0.2: ragged += 1
    cons = np.array(cons)
    ident = np.array([(B[i][G[i]] == cons[G[i]]).mean() for i in range(n)])
    med = float(np.median(ident)); mad = float(np.median(np.abs(ident - med)))
    outl = float((ident < med - max(3 * mad, 0.03)).mean())
    # variable columns: minority state (incl. gap) >= 15 % of the members
    feats = []
    for j in range(B.shape[1]):
        col = B[:, j]; c = collections.Counter(col); top = c.most_common(1)[0][1]
        for s, k in c.items():
            if k < c.most_common(1)[0][1] and k / n >= 0.15: feats.append(col == s)
    f2 = best = 0.0
    if len(feats) >= 2 and n >= 8:
        X = np.array(feats, dtype=float).T                        # members x features
        d = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1); i, j = np.unravel_index(d.argmax(), d.shape)
        cent = np.array([X[i], X[j]])
        for _ in range(15):
            lab = ((X[:, None, :] - cent[None]) ** 2).sum(-1).argmin(1)
            if len(set(lab)) < 2: break
            cent = np.array([X[lab == k].mean(0) for k in (0, 1)])
        if len(set(lab)) == 2:
            f2 = min((lab == 0).mean(), (lab == 1).mean())
            best = float(np.abs(X[lab == 0].mean(0) - X[lab == 1].mean(0)).max())
    typ = 'T2' if (f2 >= 0.2 and best >= 0.7) else ('TJ' if outl >= 0.10 else 'T1')
    return dict(n=n, id_med=med, agree=float(np.mean(agree)), ragged=ragged / B.shape[1], f2=f2, best=best, outl=outl, type=typ)

def main():
    work, prefix = sys.argv[1:3]
    truth = dict(l.split() for l in open(sys.argv[3])) if len(sys.argv) > 3 else None
    mem = collections.defaultdict(list)
    for l in open(os.path.join(work, 'chunks.tsv')):
        i, c, s = l.rstrip('\n').split('\t'); mem[c].append(i)
    print('chunk\tn\tid_med\tagree\tragged\tf2\tbest\toutl\ttype' + ('\ttruth\ttop_share' if truth else ''))
    for p in sorted(glob.glob(os.path.join(work, prefix + '_*.fasta.aln'))):
        c = os.path.basename(p)[:-len('.fasta.aln')]
        names, A = read_aln(p)
        if len(names) < 8: continue
        s = stats(A)
        if not s: continue
        row = '%s\t%d\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%s' % (c, s['n'], s['id_med'], s['agree'], s['ragged'], s['f2'], s['best'], s['outl'], s['type'])
        if truth:
            cnt = collections.Counter(truth[i] for i in mem[c]).most_common(); tot = sum(v for _, v in cnt)
            share = cnt[0][1] / tot
            tc = 'pure' if share >= 0.9 else ('two' if len(cnt) > 1 and cnt[1][1] / tot >= 0.2 and (cnt[0][1] + cnt[1][1]) / tot >= 0.9 else 'other')
            row += '\t%s\t%.2f' % (tc, share)
        print(row)

if __name__ == '__main__':
    main()
