# content-based tie-breaking for kmer_order.c (experiment): clusters are ranked by the smallest FNV hash of their sequences' content
import sys
t = open(sys.argv[1]).read()
def R(a, b):
    global t
    assert t.count(a) == 1, a[:50]; t = t.replace(a, b)
R("static float *CD;", "static uint64_t *KEY;\nstatic uint64_t fnvseq(const char *s, size_t n){uint64_t h=1469598103934665603ULL;for(size_t i=0;i<n;i++){h^=(unsigned char)s[i];h*=1099511628211ULL;}return h;}\nstatic int pair_less(size_t i1, long j1, size_t i2, long j2){uint64_t a1=KEY[i1]<KEY[j1]?KEY[i1]:KEY[j1], b1=KEY[i1]<KEY[j1]?KEY[j1]:KEY[i1], a2=KEY[i2]<KEY[j2]?KEY[i2]:KEY[j2], b2=KEY[i2]<KEY[j2]?KEY[j2]:KEY[i2];return a1<a2||(a1==a2&&b1<b2);}\nstatic float *CD;")
R("        if (ACTIVE[j] && row[j] < best) { best = row[j]; bj = (long)j; }", "        if (ACTIVE[j] && (row[j] < best || (row[j] == best && bj >= 0 && KEY[j] < KEY[bj]))) { best = row[j]; bj = (long)j; }")
R("    HH = xmalloc(n * sizeof(uint64_t));", "    KEY = xmalloc(n * sizeof(uint64_t)); for (size_t i = 0; i < n; i++) KEY[i] = fnvseq(RECS[i].seq, RECS[i].len);\n    HH = xmalloc(n * sizeof(uint64_t));")
R("if (ACTIVE[i] && NNJ[i] >= 0 && NND[i] < bd) { bd = NND[i]; ci = (long)i; }", "if (ACTIVE[i] && NNJ[i] >= 0 && (NND[i] < bd || (NND[i] == bd && ci >= 0 && pair_less(i, NNJ[i], (size_t)ci, NNJ[ci])))) { bd = NND[i]; ci = (long)i; }")
R("if (v < NND[i] || (v == NND[i] && ci < NNJ[i]))", "if (v < NND[i] || (v == NND[i] && KEY[ci] < KEY[NNJ[i]]))")
# arrangement: the cluster with the smaller key goes first; tied end-pair options by the keys of the end sequences
old_start = t.index("        int *A = cl[ci], *B = cl[cj];")
old_end = t.index("        cl[ci] = A; clen[ci] = sI + sJ; cl[cj] = NULL; clen[cj] = 0;")
new = """        int *A = cl[ci], *B = cl[cj];
        size_t sI = clen[ci], sJ = clen[cj];
        if (KEY[cj] < KEY[ci]) { int *tp = A; A = B; B = tp; size_t ts = sI; sI = sJ; sJ = ts; }
        float opts[4] = {dist_pair((size_t)A[sI - 1], (size_t)B[0]), dist_pair((size_t)A[sI - 1], (size_t)B[sJ - 1]),
                         dist_pair((size_t)A[0], (size_t)B[0]), dist_pair((size_t)A[0], (size_t)B[sJ - 1])};
        int best = 0;
        for (int t = 1; t < 4; t++) {
            if (opts[t] < opts[best]) best = t;
            else if (opts[t] == opts[best]) {
                size_t ea = (t >= 2) ? (size_t)A[0] : (size_t)A[sI - 1], eb = (t == 1 || t == 3) ? (size_t)B[sJ - 1] : (size_t)B[0];
                size_t ba = (best >= 2) ? (size_t)A[0] : (size_t)A[sI - 1], bb = (best == 1 || best == 3) ? (size_t)B[sJ - 1] : (size_t)B[0];
                if (pair_less(ea, (long)eb, ba, (long)bb)) best = t;
            }
        }
        if (best >= 2) for (size_t x = 0, y = sI - 1; x < y; x++, y--) { int t = A[x]; A[x] = A[y]; A[y] = t; }
        if (best == 1 || best == 3) for (size_t x = 0, y = sJ - 1; x < y; x++, y--) { int t = B[x]; B[x] = B[y]; B[y] = t; }
        A = realloc(A, (sI + sJ) * sizeof(int));
        if (!A) die("out of memory");
        memcpy(A + sI, B, sJ * sizeof(int));
        free(B);
        if (KEY[cj] < KEY[ci]) KEY[ci] = KEY[cj];
        { size_t sa = clen[ci], sb = clen[cj]; sI = sa; sJ = sb; }
"""
t = t[:old_start] + new + t[old_end:]
open(sys.argv[2], 'w').write(t)
