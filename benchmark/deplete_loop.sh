#!/usr/bin/env bash
# Deplete-and-resample loop for rare subfamilies, on simulated data with known truth.
#
#   1. simulate 100,000 copies: 7 abundant subfamilies + 1 rare one (100 copies, 0.1 %)
#   2. round 1: sample SAMPLE copies, SubFam -P, assign ALL copies to the round's consensuses
#   3. deplete: a copy is "explained" when its identity to its best consensus is in the upper
#      part of that consensus's identity distribution (>= the given percentile of the copies it
#      attracted); the rest (unassigned, or further from the consensus than its typical copy)
#      stays in the residual. Rare sister subfamilies sit a few % below the sister's copies, so
#      only the close half (or quarter) of each consensus's copies should be removed per round
#   4. round 2: sample again from the residual, SubFam -P; does the rare subfamily appear now?
#
#   deplete_loop.sh [outdir] [threads] [sample] [percentile]
# percentile (default 50): a copy is explained when its identity to its best consensus is at or
# above this percentile of that consensus's copies; lower keeps less, higher keeps more.
# Requires: mafft, vsearch, python3 + numpy.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-deplete_out}; T=${2:-$(nproc)}; SAMPLE=${3:-30000}; PCT=${4:-50}
mkdir -p "$OUT"; cd "$OUT"

python3 "$HERE/simulate.py" all --divergence 0.08 --sizes 30000,20000,15000,12500,10000,7500,4900,100 --seed 11
echo "truth: $(cut -f2 all.truth.tsv | sort | uniq -c | sort -k1,1nr | awk '{printf "%s=%s ", $2, $1}')"

evaluate() {  # evaluate <label> <SubFam outdir> <prefix> <sampled.fasta>
    grep '^>' "$4" | sed 's/^>//' | awk 'NR==FNR{s[$1]=1; next} ($1 in s)' - all.truth.tsv > "$2/truth.tsv"
    python3 "$HERE/evaluate.py" "$1" "$2/$3.cons.fasta" "$2/$3.chunks.tsv" "$2/truth.tsv" all.masters.fasta
}
sample() {  # sample <in.fasta> <out.fasta> <n> <seed>
    python3 - "$@" <<'PY'
import random, sys
src, dst, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
recs = open(src).read().split(">")[1:]
random.Random(seed).shuffle(recs)
open(dst, "w").write("".join(">" + r for r in recs[:n]))
PY
}

echo "== round 1"
sample all.fasta r1.fasta "$SAMPLE" 1
echo "rare copies in the sample: $(grep '^>' r1.fasta | sed 's/^>//' | awk 'NR==FNR{s[$1]=1; next} ($1 in s) && $2=="SF8"' - all.truth.tsv | wc -l)"
"$HERE/../SubFam.sh" -P -n 50 -t "$T" -o round1 -x r1 r1.fasta > round1.log
evaluate "round 1 (sample of $SAMPLE)" round1 r1 r1.fasta

echo "== assign all copies to the round-1 consensuses and deplete"
vsearch --usearch_global all.fasta --db round1/r1.cons.fasta --id 0.5 --iddef 2 --strand plus \
    --maxaccepts 4 --maxrejects 32 --threads "$T" --quiet --userout hits1.tsv --userfields query+target+id
PCT=$PCT python3 - <<'PY'
import collections, os, numpy as np
pct = float(os.environ["PCT"])
hits = {}
for line in open("hits1.tsv"):
    q, t, i = line.split("\t"); hits[q] = (t, float(i))
by_cons = collections.defaultdict(list)
for q, (t, i) in hits.items(): by_cons[t].append(i)
floor = {t: np.percentile(v, pct) for t, v in by_cons.items()}
ids = [l[1:].split()[0] for l in open("all.fasta") if l.startswith(">")]
truth = dict(l.split() for l in open("all.truth.tsv"))
explained = {q for q, (t, i) in hits.items() if i >= floor[t]}
residual = [q for q in ids if q not in explained]
with open("residual.ids", "w") as fh: fh.write("\n".join(residual) + "\n")
c_all = collections.Counter(truth[q] for q in ids); c_res = collections.Counter(truth[q] for q in residual)
print(f"explained {len(explained)} copies, residual {len(residual)}")
print("residual by subfamily: " + " ".join(f"{k}={c_res[k]}/{c_all[k]}" for k in sorted(c_all)))
PY
awk 'NR==FNR{k[$1]=1; next} /^>/{p=(substr($1,2) in k)} p' residual.ids all.fasta > residual.fasta

echo "== round 2 (from the residual)"
n=$(grep -c '^>' residual.fasta); [ "$n" -le "$SAMPLE" ] && cp residual.fasta r2.fasta || sample residual.fasta r2.fasta "$SAMPLE" 2
echo "round-2 input: $(grep -c '^>' r2.fasta) copies"
"$HERE/../SubFam.sh" -P -n 50 -t "$T" -o round2 -x r2 r2.fasta > round2.log
evaluate "round 2 (residual)" round2 r2 r2.fasta
echo "== both rounds together"
mkdir -p both
python3 - <<'PY'                       # a copy sampled in both rounds keeps its round-2 chunk
chunk = {}
for f in ("round1/r1.chunks.tsv", "round2/r2.chunks.tsv"):
    for line in open(f):
        q, c, _ = line.rstrip("\n").split("\t"); chunk[q] = c
used = set(chunk.values())
with open("both/b.chunks.tsv", "w") as fh:
    for q, c in chunk.items(): fh.write(f"{q}\t{c}\t+\n")
with open("both/b.cons.fasta", "w") as out:
    for f in ("round1/r1.cons.fasta", "round2/r2.cons.fasta"):
        keep = False
        for line in open(f):
            if line.startswith(">"): keep = line[1:].split()[0] in used
            if keep: out.write(line)
PY
cat r1.fasta r2.fasta | awk '/^>/{if(seen[$1]++) skip=1; else skip=0} !skip' > both.fasta
evaluate "rounds 1+2" both b both.fasta
