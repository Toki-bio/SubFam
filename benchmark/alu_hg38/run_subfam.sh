#!/usr/bin/env bash
# Run SubFam on a sample of hg38 Alu copies, then score consensus recovery against reference sets.
# usage: run_subfam.sh ALU_COPIES.fa SAMPLE_N SEED OUTDIR REF1.fa [REF2.fa ...]
# Needs: SubFam.sh (this repo), seqkit, mafft, python3 with numpy and biopython.
set -euo pipefail
IN=$1; N=$2; SEED=$3; OUT=$4; shift 4
HERE=$(cd "$(dirname "$0")" && pwd); SUBFAM="$HERE/../../SubFam.sh"
mkdir -p "$OUT"
seqkit sample -n "$N" -s "$SEED" "$IN" > "$OUT/sample.fa"
for n in 20 50; do                                   # chunk sizes fixed in advance (README)
    "$SUBFAM" -n "$n" -o "$OUT/n$n" -x alu "$OUT/sample.fa" > "$OUT/n$n.log" 2>&1
    for ref in "$@"; do
        python3 "$HERE/recover.py" "$ref" "$OUT/n$n/alu.cons.fasta" "$OUT/n$n.$(basename "${ref%.*}")" | sed "s|^|n=$n $(basename "$ref"): |"
    done
done
