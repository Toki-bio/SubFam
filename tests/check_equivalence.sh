#!/usr/bin/env bash
# Check that SubFam's embedded programs reproduce the tools they replace, byte for byte:
#   consensus (awk)  vs  EMBOSS cons        random alignments x plurality x setcase
#   MSF writer (awk) vs  EMBOSS seqret      random alignments, ragged lengths, long names
#   k-mer order (py) vs  ViewAlign kmer-tree.js (optional; needs node)
#
#   tests/check_equivalence.sh [path/to/MSA-viewer/kmer-tree.js]
# Requires: EMBOSS (cons, seqret), python3 with numpy; node for the third check.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SUBFAM=$HERE/../SubFam.sh
KT=${1:-}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
cd "$W"

extract() { sed -n "/^$1=\$(cat <<'$2'\$/,/^$2\$/p" "$SUBFAM" | sed '1d;$d'; }
extract CONS_AWK AWK > cons.awk
extract MSF_AWK AWK > msf.awk
extract KMER_ORDER_PY PY > order.py

python3 - <<'PY'
import random
r = random.Random(7)
alphabets = ["ACGT"] * 6 + ["ACGTN-"] * 3 + ["ACGTRYKMSWBDHVN-.X~", "acgtn-", "A-", "-"]
for t in range(200):                                       # consensus test alignments
    n = r.choice([2, 3, 4, 5, 7, 10, 20, 33, 50, 61]); L = r.randint(1, 80); a = r.choice(alphabets)
    base = [r.choice(a) for _ in range(L)]
    with open(f"c{t:03d}.fa", "w") as o:
        for i in range(n):
            keep = r.choice([.5, .8, .95])
            o.write(f">s{i}\n" + "".join(c if r.random() < keep else r.choice(a) for c in base) + "\n")
for t in range(60):                                        # MSF test alignments
    n = r.randint(1, 30); L = r.randint(1, 260)
    with open(f"m{t:02d}.fa", "w") as o:
        for i in range(n):
            ln = L if r.random() < 0.8 else r.randint(1, L)
            s = "".join(r.choice("acgtACGT-.~") if r.random() < 0.7 else "-" for _ in range(ln))
            nm = f"s{i}" if r.random() < 0.7 else f"a_much_longer_name_{i} with description"
            o.write(f">{nm}\n{s[:len(s) // 2]}\n{s[len(s) // 2:]}\n")
for L in (1, 2, 3, 49, 50, 51, 52, 53, 99, 100, 101, 102, 103, 151, 152, 153, 1001, 1002, 1003, 1004):
    with open(f"len{L:04d}.fa", "w") as o:
        o.write(f">a\n{'a' * L}\n>b\n{'-' * (L // 2)}{'c' * (L - L // 2)}\n")
PY

tot=0; bad=0
for f in c*.fa; do
    n=$(grep -c '>' "$f")
    for p in 1 2 $(( (n * 36 + 99) / 100 )) $((n / 2)) 7.5 "$n"; do
        for sc in default 3 0; do
            if [ $sc = default ]; then
                e=$(cons -sequence "$f" -outseq stdout -plurality "$p" -name X -auto 2>/dev/null)
                a=$(awk -v plur="$p" -v name=X -f cons.awk "$f")
            else
                e=$(cons -sequence "$f" -outseq stdout -plurality "$p" -setcase $sc -name X -auto 2>/dev/null)
                a=$(awk -v plur="$p" -v setcase=$sc -v name=X -f cons.awk "$f")
            fi
            tot=$((tot + 1)); [ "$e" == "$a" ] || { bad=$((bad + 1)); echo "  cons differs: $f plurality=$p setcase=$sc"; }
        done
    done
done
echo "consensus vs EMBOSS cons:  $tot runs, $bad different"
fails=$bad

tot=0; bad=0; d=$(date +%d/%m/%y)
for f in m*.fa len*.fa; do
    seqret -sequence "$f" -outseq out.msf -osformat2 msf -auto 2>/dev/null
    awk -v file=out.msf -v date="$d" -f msf.awk "$f" > mine.msf
    tot=$((tot + 1)); cmp -s out.msf mine.msf || { bad=$((bad + 1)); echo "  msf differs: $f"; }
done
echo "MSF vs EMBOSS seqret:      $tot files, $bad different"
fails=$((fails + bad))

if [ -n "$KT" ] && command -v node >/dev/null 2>&1; then
    python3 "$HERE/../benchmark/simulate.py" sim --divergence 0.10 --sizes 120,100,90,80,70,60,50,40 >/dev/null
    python3 - <<'PY'                                         # same copies, every third one reversed
seqs = open("sim.fasta").read().split(">")[1:]
comp = str.maketrans("ACGT", "TGCA")
with open("mixed.fasta", "w") as o:
    for i, b in enumerate(seqs):
        h, s = b.split("\n", 1); s = s.replace("\n", "")
        o.write(f">{h}\n{s.translate(comp)[::-1] if i % 3 == 0 else s}\n")
PY
    cat > ref.js <<'JS'
const fs = require('fs'), KT = require(process.argv[2]);
const [f, k, can] = process.argv.slice(3);
const seqs = fs.readFileSync(f, 'utf8').split('>').slice(1)
    .map(b => { const l = b.split('\n'); return { name: l[0].trim(), seq: l.slice(1).join('').trim() }; });
const t = KT.guideTree(seqs, +k, { canonical: can === '1' });
console.log(t.order.map(i => seqs[i].name).join('\n'));
JS
    tot=0; bad=0
    for f in sim.fasta mixed.fasta; do
        for k in 3 4 6 8 12; do
            for c in 0 1; do
                node ref.js "$(cd "$(dirname "$KT")" && pwd)/$(basename "$KT")" "$f" $k $c > js.txt
                python3 order.py "$f" py.fa $k $c 0; grep '>' py.fa | sed 's/^>//' > py.txt
                tot=$((tot + 1)); cmp -s js.txt py.txt || { bad=$((bad + 1)); echo "  order differs: $f k=$k canonical=$c"; }
            done
        done
    done
    echo "k-mer order vs kmer-tree.js: $tot runs, $bad different"
    fails=$((fails + bad))
else
    echo "k-mer order vs kmer-tree.js: skipped (give the path to kmer-tree.js; needs node)"
fi
[ "$fails" -eq 0 ] && echo "ALL EQUIVALENT" || { echo "DIFFERENCES FOUND"; exit 1; }
