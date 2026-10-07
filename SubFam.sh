#!/usr/bin/env bash
# SubFam: reduce a large set of related sequences (TE copies, satellite monomers,
# multigene family members) to an alignment of consensus sequences.
#
#   1. order all sequences along the MAFFT guide tree (no alignment: --retree 0 --reorder)
#   2. cut the ordered list into consecutive chunks of N neighbours
#   3. align each chunk and call its plurality consensus (EMBOSS cons)
#   4. align the consensus sequences (MAFFT L-INS-i)
#
# Each consensus averages N closely related copies: private mutations of single
# copies cancel out, variants shared by a group survive. No identity threshold
# is needed, so it works equally for young and very old (highly diverged) families.
#
# Original idea by Toki, improved by Nikita (parallel chunk alignment).

set -euo pipefail

VERSION=1.0.0
N=50            # sequences per chunk
PLURALITY=0.36  # fraction of chunk sequences that must agree for a consensus base
THREADS=$(nproc)
OUTDIR=.
PREFIX=
ADJUST=
PARTTREE=
KEEP_N=
KEEP_TMP=

usage() {
    cat <<EOF
SubFam $VERSION - reduce related sequences to an alignment of chunk consensuses

Usage: $(basename "$0") [options] <sequences.fasta>

Options:
  -n INT    sequences per chunk (default: $N)
  -p FLOAT  plurality: fraction of a chunk that must agree on a base (default: $PLURALITY)
  -t INT    threads (default: all, $THREADS)
  -o DIR    output directory (default: current directory)
  -x STR    output prefix (default: input file name without extension)
  -r        sequences may be on both strands: let MAFFT adjust direction
  -P        use MAFFT PartTree for the ordering step (recommended for >10,000 sequences)
  -a        keep no-consensus positions as N (default: drop them from the consensus)
  -k        keep intermediate files (chunks, chunk alignments)
  -v        print version
  -h        show this help

Outputs (in DIR):
  PREFIX.cons.fasta   unaligned consensus sequences, one per chunk
  PREFIX.aln.fasta    aligned consensus sequences (FASTA)
  PREFIX.msf          the same alignment in MSF format
  PREFIX.chunks.tsv   input id, consensus it went into, strand (- = reverse-complemented
                      by -r, relative to the first sequence)

Requires: mafft, cons and seqret (EMBOSS), awk.
EOF
}

die() { echo "SubFam: $*" >&2; exit 1; }
elapsed() { local s=$(( $(date +%s) - $1 )); echo "  done in $((s / 3600))h $(((s / 60) % 60))m $((s % 60))s"; }

while getopts ":n:p:t:o:x:rPakvh" opt; do
    case $opt in
        n) N=$OPTARG ;;
        p) PLURALITY=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        x) PREFIX=$OPTARG ;;
        r) ADJUST=--adjustdirection ;;
        P) PARTTREE=--parttree ;;
        a) KEEP_N=1 ;;
        k) KEEP_TMP=1 ;;
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
for tool in mafft cons seqret awk; do
    command -v "$tool" >/dev/null 2>&1 || die "'$tool' not found in PATH"
done
[[ $N =~ ^[0-9]+$ ]] && [ "$N" -ge 2 ] || die "-n must be an integer >= 2"
[[ $THREADS =~ ^[0-9]+$ ]] && [ "$THREADS" -ge 1 ] || die "-t must be a positive integer"
awk -v p="$PLURALITY" 'BEGIN { exit !(p > 0 && p <= 1) }' || die "-p must be in (0, 1]"

PREFIX=${PREFIX:-$(basename "${IN%.*}")}
NAME=${PREFIX//[[:space:]]/_}   # FASTA ids cannot contain spaces
mkdir -p "$OUTDIR"
WORK=$(mktemp -d "$OUTDIR/.subfam_${PREFIX}_XXXXXX")
[ -n "$KEEP_TMP" ] || trap 'rm -rf "$WORK"' EXIT

T0=$(date +%s)
TOTAL=$(grep -c '^>' "$IN" || true)
[ "$TOTAL" -ge 1 ] || die "no FASTA records in '$IN'"
NCHUNK=$(( TOTAL / N )); [ "$NCHUNK" -ge 1 ] || NCHUNK=1
echo "SubFam $VERSION: $TOTAL sequences -> $NCHUNK chunks of ~$N"

echo "Ordering sequences along the guide tree"
T1=$(date +%s)
mafft --thread "$THREADS" --nuc --quiet --retree 0 --reorder $PARTTREE $ADJUST "$IN" > "$WORK/ordered.fasta"
elapsed "$T1"

# Split into NCHUNK chunks of N consecutive sequences; the remainder (< N)
# joins the last chunk so that no sequence is lost.
awk -v n="$N" -v nc="$NCHUNK" -v dir="$WORK" -v pre="$NAME" -v map="$WORK/chunks.tsv" '
    BEGIN { fmt = "%s_%0" length(nc "") "d" }
    /^>/ {
        c = int(i++ / n); if (c >= nc) c = nc - 1
        name = sprintf(fmt, pre, c + 1)
        out = dir "/" name ".fasta"
        id = substr($1, 2); sub(/^_R_/, "", id)   # strip MAFFT reverse-complement tag
        print id "\t" name "\t" (substr($1, 2, 3) == "_R_" ? "-" : "+") > map
    }
    { print > out }
' "$WORK/ordered.fasta"
rm -f "$WORK/ordered.fasta"

echo "Aligning chunks and calling consensus sequences"
T1=$(date +%s)
export PLURALITY KEEP_N
find "$WORK" -maxdepth 1 -name "${NAME}_*.fasta" -print0 | sort -z |
    xargs -0 -P "$THREADS" -I {} sh -c '
        set -e; f=$1; name=$(basename "$f" .fasta)
        n=$(grep -c "^>" "$f")
        plur=$(awk -v n="$n" -v p="$PLURALITY" "BEGIN { x = n * p; print (x == int(x)) ? x : int(x) + 1 }")
        mafft --thread 1 --nuc --quiet "$f" > "$f.aln"
        cons -sequence "$f.aln" -outseq stdout -plurality "$plur" -name "$name" -auto |
            awk -v keep="$KEEP_N" "!/^>/ && !keep { gsub(/[Nn]/, \"\") } 1" > "$f.cons"
    ' _ {}
find "$WORK" -maxdepth 1 -name "${NAME}_*.fasta.cons" -print0 | sort -z | xargs -0 cat > "$OUTDIR/$PREFIX.cons.fasta"
cp "$WORK/chunks.tsv" "$OUTDIR/$PREFIX.chunks.tsv"
elapsed "$T1"

echo "Aligning consensus sequences"
T1=$(date +%s)
if [ "$NCHUNK" -ge 2 ]; then
    mafft --thread "$THREADS" --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet \
        "$OUTDIR/$PREFIX.cons.fasta" > "$OUTDIR/$PREFIX.aln.fasta"
else
    cp "$OUTDIR/$PREFIX.cons.fasta" "$OUTDIR/$PREFIX.aln.fasta"
fi
seqret -sequence "$OUTDIR/$PREFIX.aln.fasta" -outseq "$OUTDIR/$PREFIX.msf" -osformat2 msf -auto
elapsed "$T1"

[ -z "$KEEP_TMP" ] || echo "Intermediate files kept in $WORK"
echo "Job completed: $OUTDIR/$PREFIX.{cons.fasta,aln.fasta,msf,chunks.tsv}"
elapsed "$T0"
