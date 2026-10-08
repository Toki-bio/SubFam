#!/usr/bin/env bash
# Robustness of SubFam.sh on awkward but legal inputs.
#   tests/test_robust.sh [WORKDIR]
# For every case SubFam must either (a) exit 0 and keep every input record (each id in chunks.tsv exactly
# once, equal-length aligned rows, an MSF file) or (b) exit non-zero with a message on stderr. It must
# never exit 0 with sequences lost, an empty alignment, or a hang (each run has a 300 s limit).
# Needs: bash, awk, mafft, timeout, and for the default ordering a C compiler or python>=3.6+numpy.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBFAM=${SUBFAM:-$HERE/../SubFam.sh}
W=${1:-$(mktemp -d)}
mkdir -p "$W" && cd "$W" || exit 1
pass=0; fail=0
rnd() { awk -v n="$1" -v s="$2" 'BEGIN{srand(s); b="ACGT"; for(i=0;i<n;i++) printf "%s", substr(b,int(rand()*4)+1,1); print ""}'; }
mutate() { awk -v r="$2" -v s="$3" 'BEGIN{srand(s); b="ACGT"} {n=length($0); o=""; for(i=1;i<=n;i++){c=substr($0,i,1); if(rand()<r) c=substr(b,int(rand()*4)+1,1); o=o c} print o}' <<<"$1"; }
revcomp() { echo "$1" | rev | tr 'ACGTacgt' 'TGCAtgca'; }
# a family: BASE plus a second lineage with a few changes, copies with 5 % private change
BASE=$(rnd 220 1); ALT=$(mutate "$BASE" 0.03 2)
fam() {   # fam N prefix -> FASTA on stdout
    local n=$1 p=$2 i s
    for ((i = 0; i < n; i++)); do
        if ((i % 3)); then s=$(mutate "$BASE" 0.05 $((i + 100))); else s=$(mutate "$ALT" 0.05 $((i + 100))); fi
        printf '>%s%d\n%s\n' "$p" "$i" "$s"
    done
}

run_case() {   # name, expected-ids-file or "-", args...  (input is in-<name>.fa)
    local name=$1; shift
    local in=in-$name.fa out=out-$name
    rm -rf "$out"
    timeout 300 bash "$SUBFAM" "$@" -o "$out" "$in" > "log-$name.txt" 2> "err-$name.txt"
    local rc=$?
    local ok=1 why=""
    if [ $rc -eq 124 ]; then ok=0; why="timeout"
    elif [ $rc -ne 0 ]; then
        if [ -s "err-$name.txt" ]; then why="refused with a message (rc=$rc): $(head -c 120 "err-$name.txt" | tr '\n' ' ')"; else ok=0; why="rc=$rc and no message"; fi
    else
        local pre; pre=$(ls "$out"/*.chunks.tsv 2>/dev/null | head -1); pre=${pre%.chunks.tsv}
        local nin; nin=$(grep -c '^>' "$in")
        if [ -z "$pre" ]; then ok=0; why="no chunks.tsv"
        else
            local nout nuniq
            nout=$(wc -l < "$pre.chunks.tsv"); nuniq=$(cut -f1 "$pre.chunks.tsv" | sort -u | wc -l)
            [ "$nout" -eq "$nin" ] || { ok=0; why="chunks.tsv has $nout lines for $nin records"; }
            if [ $ok -eq 1 ]; then
                local widths; widths=$(awk '/^>/{if(s!="")print length(s); s=""; next}{s=s $0} END{if(s!="")print length(s)}' "$pre.aln.fasta" | sort -u | wc -l)
                [ "$widths" -le 1 ] || { ok=0; why="aligned rows differ in length"; }
                [ -s "$pre.aln.fasta" ] && [ -s "$pre.msf" ] || { ok=0; why="empty alignment or MSF"; }
            fi
            [ $ok -eq 1 ] && why="ok, $nout records in $(cut -f2 "$pre.chunks.tsv" | sort -u | wc -l) chunk(s)"
        fi
    fi
    if [ $ok -eq 1 ]; then pass=$((pass + 1)); printf 'PASS  %-22s %s\n' "$name" "$why"
    else fail=$((fail + 1)); printf 'FAIL  %-22s %s\n' "$name" "$why"; fi
}

fam 1 a > in-one.fa;                  run_case one -n 50
fam 2 a > in-two.fa;                  run_case two -n 50
fam 10 a > in-fewer_than_chunk.fa;    run_case fewer_than_chunk -n 50
fam 40 a > in-n2.fa;                  run_case n2 -n 2
fam 40 a > in-n_gt_count.fa;          run_case n_gt_count -n 1000
fam 60 a > in-plain60.fa;             run_case plain60 -n 20
{ fam 20 a; fam 20 a; } > in-duplicate_ids.fa;    run_case duplicate_ids -n 20
fam 30 a | sed 's/^>a\([0-9]*\)$/>a\1 some description with spaces|and|pipes/' > in-spaces_pipes.fa; run_case spaces_pipes -n 10
fam 30 a | awk 'NR%2==0 && NR%10==0 {print tolower($0); next} NR%2==0 && NR%14==0 {print substr($0,1,50) "NNNNNRYKM" substr($0,60); next} {print}' > in-lower_N_iupac.fa; run_case lower_N_iupac -n 10
{ fam 20 a; printf '>empty\n\n'; fam 5 b; } > in-empty_record.fa; run_case empty_record -n 10
fam 30 a | sed 's/$/\r/' > in-crlf.fa; run_case crlf -n 10
fam 30 a | fold -w 60 | awk 'NR%7==0{print ""} {print}' | head -c -1 > in-wrapped_blank_nonl.fa; run_case wrapped_blank_nonl -n 10
awk -v s="$BASE" 'BEGIN{for(i=0;i<120;i++) printf ">i%d\n%s\n", i, s}' > in-all_identical.fa; run_case all_identical -n 20
fam 30 a | awk '/^>/{print; next} {print substr($0,1,40) "---" substr($0,41)}' > in-gapped_input.fa; run_case gapped_input -n 10
fam 30 a | awk 'BEGIN{srand(3)} /^>/{h=$0; next} {print h; print substr($0,1,6)}' > in-shorter_than_k.fa; run_case shorter_than_k -n 10
{ i=0; fam 30 a | while IFS= read -r h && IFS= read -r s; do i=$((i + 1)); echo "$h"; if ((i % 2)); then echo "$s"; else revcomp "$s"; fi; done; } > in-mixed_strand.fa
run_case mixed_strand_r -n 10 -r
run_case mixed_strand_no_r -n 10
cp in-plain60.fa in-protein.fa; awk '/^>/{print;next}{gsub(/A/,"L");gsub(/C/,"K");gsub(/G/,"E");gsub(/T/,"W");print}' in-plain60.fa > in-protein.fa
run_case protein_m -n 20 -m
run_case coverage -n 20 -c
mkdir -p "dir with space"; cp in-plain60.fa "dir with space/in.fa"
rm -rf "dir with space/out"; timeout 300 bash "$SUBFAM" -n 20 -o "dir with space/out" "dir with space/in.fa" > log-space.txt 2> err-space.txt \
    && [ -s "dir with space/out/in.msf" ] && { pass=$((pass + 1)); echo "PASS  dir_with_space         ok"; } || { fail=$((fail + 1)); echo "FAIL  dir_with_space         rc/outputs ($(head -c 100 err-space.txt))"; }

# determinism: same input, different thread counts, twice
fam 120 a > in-det.fa
for t in 1 8; do for rep in 1 2; do rm -rf det-$t-$rep; bash "$SUBFAM" -n 20 -t $t -o det-$t-$rep in-det.fa > /dev/null 2>&1; done; done
same() { diff -q "$1" "$2" > /dev/null; }
for f in in-det.cons.fasta in-det.chunks.tsv in-det.aln.fasta in-det.msf; do
    if same det-1-1/$f det-1-2/$f && same det-1-1/$f det-8-1/$f && same det-8-1/$f det-8-2/$f; then
        pass=$((pass + 1)); echo "PASS  determinism $f identical over 2 runs x 2 thread counts"
    else
        fail=$((fail + 1)); echo "FAIL  determinism $f differs between runs or thread counts"
    fi
done
echo "passed $pass, failed $fail (work dir: $W)"
[ "$fail" -eq 0 ]
