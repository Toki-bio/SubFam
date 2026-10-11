#!/usr/bin/env bash
# SubFam: reduce a large set of related sequences (TE copies, satellite monomers,
# multigene family members) to an alignment of consensus sequences.
#
#   1. order all sequences along a k-mer guide tree (UPGMA; neighbours end up adjacent)
#   2. cut the ordered list into consecutive chunks of N neighbours
#   3. align each chunk (MAFFT) and call its plurality consensus
#   4. align the consensus sequences (MAFFT L-INS-i)
#
# Each consensus averages N closely related copies: private mutations of single
# copies cancel out, variants shared by a group survive. No identity threshold
# is needed, so it works equally for young and very old (highly diverged) families.
#
# Original idea by Toki, improved by Nikita (parallel chunk alignment).
# The k-mer guide tree is a port of kmer-tree.js from ViewAlign (Toki-bio/MSA-viewer);
# the consensus and MSF writers reproduce EMBOSS cons and seqret output exactly.

set -euo pipefail

VERSION=1.4.1
N=50            # sequences per chunk
PLURALITY=0.36  # fraction of chunk sequences that must agree for a consensus base
K=6            # k-mer size for the ordering tree
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
OUTDIR=.
PREFIX=
ORDER=kmer
BOTH=
PARTTREE=
KEEP_N=
KEEP_TMP=
RELCOV=
MINCOV=3       # with -c: fewest covering sequences for a consensus base
EXTRACT=       # -e: chunks whose original sequences are written to one FASTA file each
LABELS=        # -g: id<TAB>group file; the composition of every chunk by group is reported

usage() {
    cat <<EOF
SubFam $VERSION - reduce related sequences to an alignment of chunk consensuses

Usage: $(basename "$0") [options] <sequences.fasta>

Options:
  -n INT    sequences per chunk (default: $N); the input is cut into round(total/N) chunks
            of equal size, so a chunk holds between about 2/3 N and 3/2 N sequences
  -p FLOAT  plurality: fraction of a chunk that must agree on a base (default: $PLURALITY)
  -k INT    k-mer size of the ordering tree, 3-12 (default: $K, best at all divergences tested)
  -t INT    threads (default: all, $THREADS)
  -o DIR    output directory (default: current directory)
  -x STR    output prefix (default: input file name without extension)
  -r        sequences may be on both strands: strand-independent k-mers, then each copy
            is oriented like its neighbour, and the whole set like its majority strand
  -m        order with the MAFFT guide tree (mafft --retree 0 --reorder) instead of k-mers
  -P        order with MAFFT PartTree (implies -m; for very large inputs, >20,000 sequences)
  -c        coverage-relative plurality, for truncated copies (fragments, 5'-truncated LINEs):
            a base needs -p of the sequences that span that column (internal gaps count,
            end gaps do not) instead of -p of the whole chunk, with at least 3 spanning
            sequences. Default: -p of the whole chunk, which cuts a consensus down to the
            region most copies cover.
  -a        keep no-consensus positions as N (default: drop them from the consensus);
            bases supported by half the chunk or fewer are written in lower case
  -e LIST   write the ORIGINAL sequences of the listed chunks to one FASTA file each, in
            DIR/PREFIX.chunk_fasta/ (as in the input: same headers, same orientation). LIST:
            all, or comma-separated chunk numbers and ranges, e.g. 3,7,22-45 (a chunk is
            named PREFIX_NN, NN = its number in the order)
  -g FILE   labels of the input sequences chosen by you (lines: id<TAB>group; a group name
            may contain spaces when the separator is a tab; e.g. a
            RepeatMasker class or a subfamily name): reports how every chunk is composed by
            group (PREFIX.composition.tsv). SubFam does not define groups.
  -K        keep intermediate files (chunks, chunk alignments)
  -v        print version
  -h        show this help

Outputs (in DIR):
  PREFIX.cons.fasta   unaligned consensus sequences, one per chunk
  PREFIX.aln.fasta    aligned consensus sequences (FASTA)
  PREFIX.msf          the same alignment in MSF format
  PREFIX.chunks.tsv   input id, consensus it went into, strand (- = reverse-complemented by -r)
  PREFIX.chunkstats.tsv  per chunk, without any labels: n, consensus length, share of
            lower-case (weakly supported) bases, columns where at least half of the chunk has
            a base (core), agree = mean share of the covering sequences that carry the most
            frequent base of a core column, id_med and id_min = median and lowest identity of
            the chunk's sequences to that base over their core columns, ragged = share of core
            columns where at least 20 % of the chunk has a gap. Low agree / id_min / high ragged
            mean a loose chunk (a mixed group, a misalignment or outlier sequences).
  PREFIX.chunk_fasta/ (with -e)   the original sequences of the chosen chunks, one file each
  PREFIX.composition.tsv (with -g)  chunk, n, labelled, majority group, purity, count per group
  PREFIX.log          version, command line and chunk count of the run

Requires: mafft >= 7, awk, and python3 with numpy >= 1.20 (not needed with -m/-P).
Input may be gzip-compressed.
EOF
}

die() { echo "SubFam: $*" >&2; exit 1; }
elapsed() { local s=$(( $(date +%s) - $1 )); echo "  done in $((s / 3600))h $(((s / 60) % 60))m $((s % 60))s"; }

ARGS=$(printf '%q ' "$@")
while getopts ":n:p:k:t:o:x:e:g:rmPcaKvh" opt; do
    case $opt in
        n) N=$OPTARG ;;
        p) PLURALITY=$OPTARG ;;
        k) K=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        x) PREFIX=$OPTARG ;;
        e) EXTRACT=$OPTARG ;;
        g) LABELS=$OPTARG ;;
        r) BOTH=1 ;;
        m) ORDER=mafft ;;
        P) ORDER=mafft; PARTTREE=--parttree ;;
        c) RELCOV=1 ;;
        a) KEEP_N=1 ;;
        K) KEEP_TMP=1 ;;
        v) echo "SubFam $VERSION"; exit 0 ;;
        h) usage; exit 0 ;;
        :) die "option -$OPTARG needs a value" ;;
        *) die "unknown option -$OPTARG (see -h)" ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 1 ] || { usage >&2; exit 1; }
IN=$1

[ -s "$IN" ] || die "input file '$IN' not found or empty"
for tool in mafft awk; do
    command -v "$tool" >/dev/null 2>&1 || die "'$tool' not found in PATH"
done
if [ "$ORDER" = kmer ]; then
    python3 -c 'import numpy; from numpy.lib.stride_tricks import sliding_window_view' 2>/dev/null ||
        die "k-mer ordering needs python3 with numpy >= 1.20 (or use -m)"
fi
[[ $N =~ ^[0-9]+$ ]] && [ "$N" -ge 2 ] || die "-n must be an integer >= 2"
[[ $K =~ ^[0-9]+$ ]] && [ "$K" -ge 3 ] && [ "$K" -le 12 ] || die "-k must be an integer from 3 to 12"
[[ $THREADS =~ ^[0-9]+$ ]] && [ "$THREADS" -ge 1 ] || die "-t must be a positive integer"
[[ $PLURALITY =~ ^[0-9]*\.?[0-9]+$ ]] && awk -v p="$PLURALITY" 'BEGIN { exit !(p > 0 && p <= 1) }' ||
    die "-p must be a number in (0, 1]"

[ -z "$LABELS" ] || [ -s "$LABELS" ] || die "label file '$LABELS' not found or empty"
[ -z "$EXTRACT" ] || [[ $EXTRACT =~ ^(all|[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*)$ ]] ||
    die "-e needs 'all' or chunk numbers and ranges such as 3,7,22-45"

BASE=$(basename "$IN")
case $BASE in *.gz) BASE=${BASE%.gz} ;; esac
PREFIX=${PREFIX:-${BASE%.*}}
[ -n "$PREFIX" ] || die "cannot derive an output prefix from '$IN'; give one with -x"
case $PREFIX in */*|*\\*) die "-x must not contain '/' or '\\'" ;; .|..) die "-x must not be '.' or '..'" ;; esac
[ -n "${PREFIX//[[:space:]]/}" ] || die "-x must not be blank"
case $OUTDIR in *\\*) die "-o must not contain '\\'" ;; esac
NAME=${PREFIX//[[:space:]]/_}   # FASTA ids cannot contain spaces
mkdir -p "$OUTDIR"
WORK=$(mktemp -d "$OUTDIR/.subfam_${NAME}_XXXXXX")
[ -n "$KEEP_TMP" ] || trap 'rm -rf "$WORK"' EXIT
export OPENBLAS_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS

if [ "$(head -c 2 "$IN" | od -An -tx1 | tr -d ' \n')" = 1f8b ]; then     # gzip magic
    gzip -dc "$IN" > "$WORK/input.fasta" || die "cannot decompress '$IN'"
    IN=$WORK/input.fasta
fi
if [ "$(head -c 3 "$IN" | od -An -tx1 | tr -d ' \n')" = efbbbf ]; then   # UTF-8 byte-order mark: would hide the first header
    tail -c +4 "$IN" > "$WORK/input_nobom.fasta"; IN=$WORK/input_nobom.fasta
fi
[ "$(head -c 1 "$IN")" = ">" ] || die "'$1' does not start with a FASTA header line"

# ---------------------------------------------------------------------------------------
# Embedded programs
# ---------------------------------------------------------------------------------------
KMER_ORDER_PY=$(cat <<'PY'
import sys
import numpy as np

# k-mer guide-tree order: a port of kmer-tree.js (ViewAlign, Toki-bio/MSA-viewer, MIT).
# Weighted Jaccard distance on k-mer counts, UPGMA with ViewAlign's tie-breaking, and at
# every merge the two clusters are flipped so that their closest ends meet.
#   args: in.fasta out.fasta K CANONICAL(0/1) ORIENT(0/1)
# Reverse-complemented sequences (ORIENT=1) get MAFFT's "_R_" prefix on their id.

CODE = np.full(256, -1, dtype=np.int64)
for ch, v in zip("ACGTUacgtu", (0, 1, 2, 3, 3, 0, 1, 2, 3, 3)):
    CODE[ord(ch)] = v


def read_fasta(path):
    names, seqs, cur = [], [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip())
                cur = []
                seqs.append(cur)
            elif cur is not None:
                cur.append(line.strip())
    return names, ["".join(s) for s in seqs]


def kmer_codes(seq, k, canonical):
    """k-mer codes over the A/C/G/T/U letters of seq (others skipped, as in kmer-tree.js)."""
    v = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    v = v[v >= 0]
    if len(v) < k:
        return np.zeros(0, dtype=np.int64)
    w = np.lib.stride_tricks.sliding_window_view(v, k)
    h = (w << (2 * np.arange(k - 1, -1, -1))).sum(1)
    if canonical:
        r = ((3 - w) << (2 * np.arange(k))).sum(1)
        h = np.where(r < h, r, h)
    return h


def distances(seqs, k, canonical):
    """float32 weighted Jaccard distance matrix: 1 - sum(min) / sum(max) of k-mer counts."""
    n = len(seqs)
    codes = [kmer_codes(s, k, canonical) for s in seqs]
    total = np.array([len(c) for c in codes], dtype=np.float64)
    # sparse count table over k-mers seen in >= 2 sequences (others never contribute to sum(min))
    rows = np.concatenate([np.full(len(c), i) for i, c in enumerate(codes)]) if n else np.zeros(0, int)
    allc = np.concatenate(codes) if n else np.zeros(0, int)
    pair = np.unique(rows * (1 << 26) + np.searchsorted(np.unique(allc), allc), return_counts=True)
    r_idx, c_idx, cnt = pair[0] >> 26, pair[0] & ((1 << 26) - 1), pair[1]
    seen = np.bincount(c_idx)
    keep = seen[c_idx] >= 2
    r_idx, c_idx, cnt = r_idx[keep], c_idx[keep], cnt[keep]
    _, c_idx = np.unique(c_idx, return_inverse=True)
    m = int(c_idx.max()) + 1 if len(c_idx) else 0
    shared = np.zeros((n, n), dtype=np.float32)          # integer counts, exact below 2^24
    # sum(min(a, b)) = sum over t >= 1 of [a >= t][b >= t]: a matrix product per count level
    block = max(1, min(m, (1 << 24) // max(n, 1)))
    for t in range(1, int(cnt.max()) + 1 if len(cnt) else 1):
        sel = cnt >= t
        rr, cc = r_idx[sel], c_idx[sel]
        if not len(cc):
            break
        _, cc = np.unique(cc, return_inverse=True)       # only the k-mers still at this count level
        mt = int(cc.max()) + 1
        for lo in range(0, mt, block):
            s = (cc >= lo) & (cc < lo + block)
            if not s.any():
                continue
            B = np.zeros((n, block), dtype=np.float32)
            B[rr[s], cc[s] - lo] = 1
            for r in range(0, n, 1024):                  # row blocks: no n x n temporary
                shared[r:r + 1024] += B[r:r + 1024] @ B.T
    d = shared                                            # converted in place, a row block at a time
    for lo in range(0, n, 256):
        sh = shared[lo:lo + 256].astype(np.float64)
        union = total[lo:lo + 256, None] + total[None, :] - sh
        with np.errstate(invalid="ignore", divide="ignore"):
            d[lo:lo + 256] = np.where(union > 0, 1 - sh / union, 1.0)
    np.fill_diagonal(d, 0)
    profiles = [np.unique(c, return_counts=True) for c in codes]
    return d, profiles, total


def pair_distance(profiles, total, a, b):
    """One entry of the distance matrix, recomputed exactly (float32, as stored)."""
    (ka, ca), (kb, cb) = profiles[a], profiles[b]
    _, ia, ib = np.intersect1d(ka, kb, assume_unique=True, return_indices=True)
    shared = float(np.minimum(ca[ia], cb[ib]).sum())
    union = total[a] + total[b] - shared
    return np.float32(1 - shared / union if union > 0 else 1.0) if a != b else np.float32(0)


def upgma_order(cd, dist):
    """Leaf order of kmer-tree.js upgma() (average linkage), same tie-breaking.
    cd: distance matrix, overwritten; dist(a, b): original distance of two sequences."""
    n = len(cd)
    if n < 2:
        return list(range(n))
    clusters = [[i] for i in range(n)]
    active = np.ones(n, dtype=bool)
    nnD = np.full(n, np.inf)
    nnJ = np.full(n, -1)

    def scan(i):
        row = np.where(active[i + 1:], cd[i, i + 1:], np.float32(np.inf))
        if len(row) and np.isfinite(row.min()):
            j = int(row.argmin())
            nnD[i], nnJ[i] = row[j], i + 1 + j
        else:
            nnD[i], nnJ[i] = np.inf, -1

    for i in range(n):
        scan(i)
    idx = np.arange(n)
    for _ in range(n - 1):
        cand = np.where(active & (nnJ >= 0), nnD, np.inf)
        ci = int(cand.argmin())
        if not np.isfinite(cand[ci]):
            break
        cj = int(nnJ[ci])
        A, B = clusters[ci], clusters[cj]
        opts = [dist(A[-1], B[0]), dist(A[-1], B[-1]), dist(A[0], B[0]), dist(A[0], B[-1])]
        best = int(np.argmin(opts))                      # first minimum, like a stable sort
        A = A[::-1] if best >= 2 else A
        B = B[::-1] if best in (1, 3) else B
        sI, sJ = len(A), len(B)
        clusters[ci], clusters[cj] = A + B, None
        active[cj] = False
        upd = active.copy()
        upd[ci] = False
        new = ((cd[ci].astype(np.float64) * sI + cd[cj].astype(np.float64) * sJ) / (sI + sJ)).astype(np.float32)
        cd[ci, upd] = new[upd]
        cd[upd, ci] = new[upd]
        scan(ci)
        others = active.copy()
        others[ci] = False
        rescan = others & ((nnJ == ci) | (nnJ == cj))
        for i in np.flatnonzero(rescan):
            scan(int(i))
        low = others & ~rescan & (idx < ci)
        v = cd[low, ci].astype(np.float64)
        win = (v < nnD[low]) | ((v == nnD[low]) & (ci < nnJ[low]))
        li = np.flatnonzero(low)[win]
        nnD[li], nnJ[li] = v[win], ci
    return next(c for i, c in enumerate(clusters) if active[i])


COMP = str.maketrans("ACGTUacgtuRYKMBVDHrykmbvdh", "TGCAAtgcaaYRMKVBHDyrmkvbhd")


def main():
    src, out, k, canonical, orient = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] == "1", sys.argv[5] == "1"
    names, seqs = read_fasta(src)
    d, profiles, total = distances(seqs, k, canonical)
    order = upgma_order(d, lambda a, b: pair_distance(profiles, total, a, b))
    flip = [False] * len(seqs)
    if orient:
        # orient each sequence like its predecessor in the order: its nearest neighbour
        prev = None
        for i in order:
            if prev is not None:
                ref = set(kmer_codes(prev, k, False).tolist())
                fwd = sum(c in ref for c in kmer_codes(seqs[i], k, False).tolist())
                rev = sum(c in ref for c in kmer_codes(seqs[i].translate(COMP)[::-1], k, False).tolist())
                flip[i] = rev > fwd
            prev = seqs[i].translate(COMP)[::-1] if flip[i] else seqs[i]
        if sum(flip) * 2 > len(flip):                    # keep the majority strand as given
            flip = [not f for f in flip]
    with open(out, "w") as fh:
        for i in order:
            if flip[i]:
                fh.write(f">_R_{names[i]}\n{seqs[i].translate(COMP)[::-1]}\n")
            else:
                fh.write(f">{names[i]}\n{seqs[i]}\n")


if __name__ == "__main__":
    main()
PY
)

CONS_AWK=$(cat <<'AWK'
# Plurality consensus of an aligned nucleotide FASTA: EMBOSS cons (embConsCalc) with its
# defaults (EDNAFULL matrix, unit weights), reimplemented. Variables: plur (plurality,
# absolute count), setcase (default: half the number of sequences), name.
# rel=1 (SubFam -c): the threshold is frac * (sequences spanning the column) instead of plur,
# needs mincov spanning sequences, and end gaps do not take part in the base choice.
BEGIN {
    sym = "A T G C S W R Y K M B V H D N U"
    nsym = split(sym, S, " ")
    split("5 -4 -4 -4 -4 1 1 -4 -4 1 -4 -1 -1 -1 -2 -4 " \
          "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5 " \
          "-4 -4 5 -4 1 -4 1 -4 1 -4 -1 -1 -4 -1 -2 -4 " \
          "-4 -4 -4 5 1 -4 -4 1 -4 1 -1 -1 -1 -4 -2 -4 " \
          "-4 -4 1 1 -1 -4 -2 -2 -2 -2 -1 -1 -3 -3 -1 -4 " \
          "1 1 -4 -4 -4 -1 -2 -2 -2 -2 -3 -3 -1 -1 -1 1 " \
          "1 -4 1 -4 -2 -2 -1 -4 -2 -2 -3 -1 -3 -1 -1 -4 " \
          "-4 1 -4 1 -2 -2 -4 -1 -2 -2 -1 -3 -1 -3 -1 1 " \
          "-4 1 1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1 1 " \
          "1 -4 -4 1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1 -4 " \
          "-4 -1 -1 -1 -1 -3 -3 -1 -1 -3 -1 -2 -2 -2 -1 -1 " \
          "-1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1 -4 " \
          "-1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1 -1 " \
          "-1 -1 -1 -4 -3 -1 -1 -3 -1 -3 -2 -2 -2 -1 -1 -1 " \
          "-2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 " \
          "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5", V, " ")       # EDNAFULL
    for (a = 1; a <= nsym; a++) {
        isym[S[a]] = 1
        for (b = 1; b <= nsym; b++) M[S[a], S[b]] = V[(a - 1) * nsym + b]
    }
}
/^>/ { n++; next }
{ gsub(/[ \t\r]/, ""); gsub(/[.~]/, "-"); gsub(/[Xx]/, "N"); seq[n] = seq[n] toupper($0) }   # as EMBOSS reads DNA
END {
    if (setcase == "") setcase = n / 2
    for (i = 1; i <= n; i++) if (length(seq[i]) > L) L = length(seq[i])
    if (rel) for (i = 1; i <= n; i++) {          # span of each sequence: first..last residue
        t = seq[i]; fi[i] = 0; la[i] = -1
        if (match(t, /[^-]/)) { fi[i] = RSTART; sub(/-+$/, "", t); la[i] = length(t) }
    }
    out = ""
    for (k = 1; k <= L; k++) {
        delete cnt; ncov = 0
        for (i = 1; i <= n; i++) {
            c = substr(seq[i], k, 1); col[i] = c
            if (c in isym) cnt[c]++
            if (rel && k >= fi[i] && k <= la[i]) ncov++
        }
        if (rel) {
            need = frac * ncov; x = int(need + 1e-9); need = (need - x > 1e-9) ? x + 1 : x   # ceil, safe for 0.07 * 100
            mc = (mincov < n) ? mincov : n
            if (ncov < mc) { out = out "N"; continue }
        }
        delete sc; delete mt
        for (s in cnt) {                       # score of one sequence with symbol s, and its +ve matches
            x = -M[s, s]; y = 0
            for (c in cnt) { x += cnt[c] * M[s, c]; if (M[s, c] > 0) y += cnt[c] }
            sc[s] = x; mt[s] = y
        }
        hi = 0; max = "unset"
        for (i = 1; i <= n; i++) {             # first sequence with the top score; a gap yields ties
            if (rel && !(col[i] in isym) && (k < fi[i] || k > la[i])) continue   # end gap: not a candidate
            v = (col[i] in isym) ? sc[col[i]] : 0
            if (max == "unset" || v > max || (v == max && col[hi] == "-")) { hi = i; max = v }
        }
        r = col[hi]; m = (r in isym) ? mt[r] : 0
        res = (m >= (rel ? need : plur)) ? r : "N"
        if (m <= setcase) res = tolower(res)
        out = out res
    }
    print ">" name
    for (i = 1; i <= length(out); i += 60) print substr(out, i, 60)
}
AWK
)

STATS_AWK=$(cat <<'AWK'
# Label-free quality of one chunk, from its alignment (FASTA, all rows the same length).
# Variables: name, conslen, lower. Prints one tab-separated line:
# chunk n conslen lower core agree id_med id_min ragged
/^>/ { n++; next }
{ gsub(/[ \t\r]/, ""); seq[n] = seq[n] toupper($0) }
END {
    L = length(seq[1]); core = 0; agsum = 0; rag = 0
    for (k = 1; k <= L; k++) {
        delete cnt; ng = 0
        for (i = 1; i <= n; i++) { c = substr(seq[i], k, 1); if (c ~ /[ACGT]/) { cnt[c]++; ng++ } }
        if (ng < n / 2) continue
        best = ""; bk = 0
        for (c in cnt) if (cnt[c] > bk || (cnt[c] == bk && c < best)) { best = c; bk = cnt[c] }
        core++; top[k] = best; agsum += bk / ng
        if ((n - ng) / n >= 0.2) rag++
    }
    for (i = 1; i <= n; i++) {
        m = 0; t = 0
        for (k = 1; k <= L; k++) if (k in top) { c = substr(seq[i], k, 1); if (c ~ /[ACGT]/) { t++; if (c == top[k]) m++ } }
        id[i] = t ? m / t : 0
    }
    for (i = 2; i <= n; i++) { v = id[i]; j = i - 1; while (j >= 1 && id[j] > v) { id[j + 1] = id[j]; j-- } id[j + 1] = v }
    med = (n % 2) ? id[(n + 1) / 2] : (id[n / 2] + id[n / 2 + 1]) / 2
    printf "%s\t%d\t%d\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", name, n, conslen, conslen ? lower / conslen : 0, core, core ? agsum / core : 0, med, id[1], core ? rag / core : 0
}
AWK
)

MSF_AWK=$(cat <<'AWK'
# Aligned FASTA -> GCG MSF, as EMBOSS seqret -osformat2 msf writes it (nucleotides).
# Variables: file (name printed in the header), date (dd/mm/yy).
function check(s,   i, c, t) {                 # GCG checksum
    t = 0; c = 0
    for (i = 1; i <= length(s); i++) { if (++c > 57) c = 1; t += c * ORD[toupper(substr(s, i, 1))] }
    return t % 10000
}
BEGIN { for (i = 32; i < 127; i++) ORD[sprintf("%c", i)] = i }
/^>/ { n++; split(substr($0, 2), w, /[ \t]/); name[n] = w[1]; next }
{ gsub(/[ \t\r]/, ""); seq[n] = seq[n] $0 }
END {
    L = 0; W = 10
    for (i = 1; i <= n; i++) { if (length(seq[i]) > L) L = length(seq[i]); if (length(name[i]) > W) W = length(name[i]) }
    tot = 0
    for (i = 1; i <= n; i++) {
        s = seq[i]
        gsub(/[ .~-]/, ".", s)                   # internal gaps '.', end gaps and padding '~'
        if (match(s, /^\.+/)) s = sprintf("%" RLENGTH "s", "") substr(s, RLENGTH + 1)
        if (match(s, /\.+$/)) s = substr(s, 1, RSTART - 1) sprintf("%" RLENGTH "s", "")
        while (length(s) < L) s = s " "
        gsub(/ /, "~", s)
        seq[i] = s; ck[i] = check(s); tot = (tot + ck[i]) % 10000
    }
    printf "!!NA_MULTIPLE_ALIGNMENT 1.0\n\n"
    printf "  %s MSF: %d Type: N %s CompCheck: %4d ..\n\n", file, L, date, tot
    for (i = 1; i <= n; i++) printf "  Name: %-" W "s Len: %d  Check: %4d Weight: 1.00\n", name[i], length(seq[i]), ck[i]
    printf "\n//\n\n"
    for (p = 1; p <= L; p += 50) {
        e = p + 49; if (e > L) e = L
        if (e == L) {
            g = e - p - length(p "")
            if (g < 0) printf "%" W "s %d %d\n", " ", p, e      # EMBOSS: unsigned width wraps
            else if (g >= length(e "")) printf "%" W "s %d %" g "d\n", " ", p, e
            else printf "           %d\n", p
        } else printf "           %-25d%25d\n", p, e
        for (i = 1; i <= n; i++) printf "%-" W "s %s\n", name[i], substr(seq[i], p, e - p + 1)
        printf "\n"
    }
}
AWK
)

# ---------------------------------------------------------------------------------------
T0=$(date +%s)
TOTAL=$(grep -c '^>' "$IN" || true)
[ "$TOTAL" -ge 1 ] || die "no FASTA records in '$IN'"
NCHUNK=$(( (TOTAL + N / 2) / N )); [ "$NCHUNK" -ge 1 ] || NCHUNK=1     # round(TOTAL / N)
if [ "$ORDER" = kmer ] && [ "$TOTAL" -gt 20000 ]; then
    echo "warning: the k-mer ordering holds a ${TOTAL}x${TOTAL} float32 matrix (~$(( 4 * TOTAL * TOTAL / 1000000000 + 1 )) GB); use -P beyond the memory of this machine" >&2
fi
echo "SubFam $VERSION: $TOTAL sequences -> $NCHUNK chunks of ~$(( TOTAL / NCHUNK ))"
{
    echo "SubFam $VERSION  $(date '+%Y-%m-%d %H:%M')"
    echo "command: $0 $ARGS"
    echo "input: $IN  sequences: $TOTAL  chunks: $NCHUNK  n: $N  p: $PLURALITY  k: $K  order: $ORDER${PARTTREE:+ parttree}${BOTH:+ both-strands}${RELCOV:+ coverage-relative}"
} > "$OUTDIR/$PREFIX.log"

T1=$(date +%s)
if [ "$ORDER" = kmer ]; then
    echo "Ordering sequences along the $K-mer guide tree"
    python3 -c "$KMER_ORDER_PY" "$IN" "$WORK/ordered.fasta" "$K" "${BOTH:-0}" "${BOTH:-0}"
else
    echo "Ordering sequences along the MAFFT guide tree"
    mafft --thread "$THREADS" --nuc --quiet --retree 0 --reorder $PARTTREE ${BOTH:+--adjustdirection} \
        "$IN" > "$WORK/ordered.fasta"
fi
elapsed "$T1"

# Split into NCHUNK consecutive chunks of equal size (within one sequence of each other).
awk -v tot="$TOTAL" -v nc="$NCHUNK" -v dir="$WORK" -v pre="$NAME" -v map="$WORK/chunks.tsv" '
    BEGIN { fmt = "%s_%0" length(nc "") "d" }
    /^>/ {
        c = int(i * nc / tot); i++
        name = sprintf(fmt, pre, c + 1)
        f = dir "/" name ".fasta"
        if (f != out) { if (out) close(out); out = f }
        id = substr($1, 2); sub(/^_R_/, "", id)   # strip the reverse-complement tag
        print id "\t" name "\t" (substr($1, 2, 3) == "_R_" ? "-" : "+") > map
    }
    { print > out }
' "$WORK/ordered.fasta"
rm -f "$WORK/ordered.fasta"

echo "Aligning chunks and calling consensus sequences"
T1=$(date +%s)
export PLURALITY KEEP_N CONS_AWK STATS_AWK RELCOV MINCOV
printf '%s\0' "$WORK/${NAME}_"*.fasta |
    xargs -0 -P "$THREADS" -I {} sh -c '
        f=$1; name=$(basename "$f" .fasta)
        n=$(grep -c "^>" "$f")
        plur=$(awk -v n="$n" -v p="$PLURALITY" "BEGIN { x = n * p; i = int(x + 1e-9); print (x - i > 1e-9) ? i + 1 : i }")
        if [ "$n" -ge 2 ]; then mafft --thread 1 --nuc --quiet "$f" > "$f.aln"; else cp "$f" "$f.aln"; fi
        [ -s "$f.aln" ] || { echo "SubFam: MAFFT failed on chunk $name" >&2; exit 1; }
        awk -v plur="$plur" -v name="$name" -v rel="$RELCOV" -v frac="$PLURALITY" -v mincov="$MINCOV" "$CONS_AWK" "$f.aln" > "$f.raw"
        [ -s "$f.raw" ] || { echo "SubFam: consensus failed on chunk $name" >&2; exit 1; }
        awk -v keep="$KEEP_N" "!/^>/ { if (keep) gsub(/n/, \"N\"); else gsub(/[Nn]/, \"\") } 1" "$f.raw" > "$f.cons"
        cl=$(grep -v "^>" "$f.cons" | tr -d "\n" | wc -c); lw=$(grep -v "^>" "$f.cons" | tr -d "\n" | tr -cd "[:lower:]" | wc -c)
        awk -v name="$name" -v conslen="$cl" -v lower="$lw" "$STATS_AWK" "$f.aln" > "$f.stats"
    ' _ {}
cat "$WORK/${NAME}_"*.fasta.cons > "$OUTDIR/$PREFIX.cons.fasta"
[ "$(grep -c '^>' "$OUTDIR/$PREFIX.cons.fasta")" -eq "$NCHUNK" ] || die "expected $NCHUNK consensus sequences, see $WORK"
cp "$WORK/chunks.tsv" "$OUTDIR/$PREFIX.chunks.tsv"
{ printf 'chunk\tn\tcons_len\tlower\tcore\tagree\tid_med\tid_min\tragged\n'; cat "$WORK/${NAME}_"*.fasta.stats; } > "$OUTDIR/$PREFIX.chunkstats.tsv"
elapsed "$T1"

if [ -n "$EXTRACT" ]; then
    # the original sequences of the chosen chunks, one FASTA file per chunk (input headers and orientation unchanged)
    XDIR="$OUTDIR/$PREFIX.chunk_fasta"; mkdir -p "$XDIR"
    awk -v sel="$EXTRACT" -v nc="$NCHUNK" -v pre="$NAME" -v dir="$XDIR" -v map="$WORK/chunks.tsv" '
        BEGIN {
            fmt = "%s_%0" length(nc "") "d"
            if (sel == "all") { for (c = 1; c <= nc; c++) want[sprintf(fmt, pre, c)] = 1 }
            else { m = split(sel, part, ","); for (q = 1; q <= m; q++) { r = split(part[q], ab, "-"); lo = ab[1] + 0; hi = (r > 1) ? ab[2] + 0 : lo
                   for (c = lo; c <= hi; c++) if (c >= 1 && c <= nc) want[sprintf(fmt, pre, c)] = 1 } }
            while ((getline line < map) > 0) { split(line, f, "\t"); if (f[2] in want) dest[f[1]] = f[2] }
        }
        /^>/ { id = substr($1, 2); out = (id in dest) ? dir "/" dest[id] ".fasta" : ""; if (out != "") nseq[dest[id]]++ }
        out != "" { print > out }
        END { for (c in want) { k++; tot += nseq[c] } printf "  wrote %d chunk files with %d sequences to %s\n", k, tot, dir
              if (!k) { print "SubFam: -e " sel " names no chunk of 1.." nc > "/dev/stderr"; exit 1 } }
    ' "$IN"
fi
if [ -n "$LABELS" ]; then
    # composition of every chunk by the groups given in the label file (SubFam only counts; it does not define groups)
    awk -F'\t' -v lab="$LABELS" '
        BEGIN { while ((getline line < lab) > 0) { sub(/\r$/, "", line); if (line ~ /^#/ || line ~ /^[ \t]*$/) continue
                    if (index(line, "\t")) split(line, a, "\t"); else split(line, a, / +/); g[a[1]] = a[2]; seen[a[2]] = 1 } ng = 0; for (x in seen) grp[++ng] = x
                for (i = 2; i <= ng; i++) { v = grp[i]; j = i - 1; while (j >= 1 && grp[j] > v) { grp[j + 1] = grp[j]; j-- } grp[j + 1] = v } }
        { n[$2]++; if (!($2 in first)) { first[$2] = 1; ord[++nch] = $2 } if ($1 in g) { cnt[$2, g[$1]]++; lb[$2]++ } }
        END {
            printf "chunk\tn\tlabelled\tmajority\tpurity"; for (i = 1; i <= ng; i++) printf "\t%s", grp[i]; printf "\n"
            for (c = 1; c <= nch; c++) { ch = ord[c]; maj = "-"; mx = 0
                for (i = 1; i <= ng; i++) if (cnt[ch, grp[i]] > mx) { mx = cnt[ch, grp[i]]; maj = grp[i] }
                printf "%s\t%d\t%d\t%s\t%s", ch, n[ch], lb[ch], maj, lb[ch] ? sprintf("%.2f", mx / lb[ch]) : "NA"
                for (i = 1; i <= ng; i++) printf "\t%d", cnt[ch, grp[i]]; printf "\n" }
        }
    ' "$OUTDIR/$PREFIX.chunks.tsv" > "$OUTDIR/$PREFIX.composition.tsv"
    echo "  chunk composition by $(head -1 "$OUTDIR/$PREFIX.composition.tsv" | awk -F'\t' '{ print NF - 5 }') groups: $OUTDIR/$PREFIX.composition.tsv"
fi

echo "Aligning consensus sequences"
T1=$(date +%s)
if [ "$NCHUNK" -ge 2 ]; then
    # L-INS-i up to 2,000 consensus sequences, MAFFT's automatic choice beyond; --threadit 0 keeps
    # the iterative refinement deterministic whatever the thread count
    MODE="--localpair --maxiterate 1000"; [ "$NCHUNK" -le 2000 ] || MODE=--auto
    mafft --thread "$THREADS" --threadit 0 $MODE --nuc --reorder --preservecase --quiet \
        "$OUTDIR/$PREFIX.cons.fasta" > "$OUTDIR/$PREFIX.aln.fasta"
else
    cp "$OUTDIR/$PREFIX.cons.fasta" "$OUTDIR/$PREFIX.aln.fasta"
fi
awk -v file="$OUTDIR/$PREFIX.msf" -v date="$(date +%d/%m/%y)" "$MSF_AWK" "$OUTDIR/$PREFIX.aln.fasta" > "$OUTDIR/$PREFIX.msf"
elapsed "$T1"

[ -z "$KEEP_TMP" ] || echo "Intermediate files kept in $WORK"
echo "Job completed: $OUTDIR/$PREFIX.{cons.fasta,aln.fasta,msf,chunks.tsv,chunkstats.tsv}"
elapsed "$T0"
