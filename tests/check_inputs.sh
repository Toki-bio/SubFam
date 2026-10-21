#!/usr/bin/env bash
# check_inputs.sh [WORKDIR] - input handling of SubFam.sh on small synthetic files (seed 1); exits 1 on any failure.
# Cases: prefix "input" (the name SINEderella uses) for plain, gzip and BOM input; CRLF; duplicate ids; empty records;
# protein; header without id; an id starting with _R_; too little free memory for the k-mer matrix; the -a
# no-consensus warning written to the log; a module named numpy.py in the current directory.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); SF="$HERE/../SubFam.sh"
W=${1:-$(mktemp -d)}; mkdir -p "$W"; cd "$W" || exit 1
python3 - <<'PY'
import random
random.seed(1)
b = "".join(random.choice("ACGT") for _ in range(200))
mut = lambda s, r: "".join(c if random.random() > r else random.choice("ACGT") for c in s)
open("input.fasta", "w").write("".join(f">s{i} some description\n{mut(b, .15)}\n" for i in range(120)))
open("dup.fa", "w").write("".join(f">x\n{mut(b, .1)}\n" for i in range(6)))
open("empty.fa", "w").write(">a\n\n" + "".join(f">s{i}\n{mut(b, .1)}\n" for i in range(5)))
open("prot.fa", "w").write("".join(f">p{i}\nMKLVQEEFFHHPPLL\n" for i in range(5)))
open("noid.fa", "w").write(">\nACGTACGT\n")
open("rr.fa", "w").write("".join(f">{'_R_' if i == 3 else ''}r{i}\n{mut(b, .1)}\n" for i in range(6)))
open("rand.fa", "w").write("".join(f">u{i}\n{''.join(random.choice('ACGT') for _ in range(200))}\n" for i in range(100)))
open("crlf.fa", "w", newline="").write("".join(f">c{i}\r\n{mut(b, .1)}\r\n" for i in range(5)))
PY
gzip -kf input.fasta; printf '\xef\xbb\xbf' | cat - input.fasta > bom.fa
fail=0
ok()  { if bash "$SF" "$@" > log.txt 2>&1; then echo "PASS  runs:  $*"; else echo "FAIL  runs:  $* ($(tail -1 log.txt))"; fail=1; fi; }
bad() { local want=$1; shift; if bash "$SF" "$@" > log.txt 2>&1; then echo "FAIL  must stop: $*"; fail=1
        elif grep -q "$want" log.txt; then echo "PASS  stops: $* ($want)"; else echo "FAIL  wrong message: $* ($(tail -1 log.txt))"; fail=1; fi; }
ok -o plain input.fasta
ok -o gz input.fasta.gz
ok -o bom -x input bom.fa
ok -o crlf crlf.fa
bad 'without sequence' -o emp empty.fa
bad 'duplicate sequence id' -o dup dup.fa
bad 'nucleotide' -o prot prot.fa
bad 'has no id' -o noid noid.fa
bad 'starts with _R_' -o rr rr.fa
SUBFAM_AVAIL_MB=1 bad 'needs ~' -o mem input.fasta
ok -a -o nfrac rand.fa
grep -q 'have no consensus (N)' nfrac/rand.log && echo "PASS  the -a no-consensus warning is in the log" \
    || { echo "FAIL  the -a no-consensus warning is not in the log"; fail=1; }
cmp -s plain/input.cons.fasta gz/input.cons.fasta && cmp -s plain/input.aln.fasta bom/input.aln.fasta \
    && echo "PASS  plain, gzip and BOM input give the same output" || { echo "FAIL  plain/gzip/BOM outputs differ"; fail=1; }
mkdir -p hij && ( cd hij && echo 'raise SystemExit("HIJACKED")' > numpy.py && bash "$SF" -o o ../input.fasta > ../log.txt 2>&1 ) \
    && echo "PASS  numpy.py in the current directory is not imported" || { echo "FAIL  numpy.py in the current directory was imported"; fail=1; }
[ $fail -eq 0 ] && echo "ALL INPUT CHECKS PASSED" || echo "INPUT CHECKS FAILED"
exit $fail
