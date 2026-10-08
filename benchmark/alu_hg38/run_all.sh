#!/usr/bin/env bash
# One-shot run of the hg38 Alu consensus-recovery protocol (see README.md; metrics are fixed there).
# usage: run_all.sh HG38_FA RMSK_TXT_GZ DFAM_ALU_FA PRICE_SUPPLEMENT_TXT OUTDIR
# Needs: bedtools, seqkit, mafft, vsearch, python3 with numpy + biopython. UNTESTED on hg38 (parts tested on small data).
set -euo pipefail
FA=$1; RMSK=$2; DFAM=$3; PRICE=$4; OUT=$5
HERE=$(cd "$(dirname "$0")" && pwd)
for t in bedtools seqkit mafft vsearch python3; do command -v $t >/dev/null || { echo "missing: $t"; exit 1; }; done
mkdir -p "$OUT"; date > "$OUT/run_date.txt"
"$HERE/extract_alu.sh" "$FA" "$RMSK" "$OUT/hg38_alu_full.fa" | tee "$OUT/extract.log"
python3 "$HERE/prep_price.py" "$PRICE" > "$OUT/price213.fa"
for seed in 1 2 3; do
  "$HERE/run_subfam.sh" "$OUT/hg38_alu_full.fa" 30000 $seed "$OUT/s$seed" "$DFAM" "$OUT/price213.fa" | tee "$OUT/s$seed.summary.txt"
  S="$OUT/s$seed/sample.fa"
  for id in 0.90 0.95 0.98; do           # baseline (a): VSEARCH centroids
    vsearch --cluster_fast "$S" --id $id --centroids "$OUT/s$seed/vs_$id.fa" --threads 8 --quiet
    for ref in "$DFAM" "$OUT/price213.fa"; do
      python3 "$HERE/recover.py" "$ref" "$OUT/s$seed/vs_$id.fa" "$OUT/s$seed/vs_$id.$(basename "${ref%.*}")" | sed "s|^|vsearch id=$id $(basename "$ref"): |" | tee -a "$OUT/s$seed.summary.txt"
    done
  done
  for n in 20 50; do                       # baseline (b): MAFFT guide-tree ordering
    "$HERE/../../SubFam.sh" -m -n $n -o "$OUT/s$seed/m$n" -x alu "$S" > "$OUT/s$seed/m$n.log" 2>&1
    for ref in "$DFAM" "$OUT/price213.fa"; do
      python3 "$HERE/recover.py" "$ref" "$OUT/s$seed/m$n/alu.cons.fasta" "$OUT/s$seed/m$n.$(basename "${ref%.*}")" | sed "s|^|subfam -m n=$n $(basename "$ref"): |" | tee -a "$OUT/s$seed.summary.txt"
    done
  done
done
echo "done; summaries: $OUT/s*.summary.txt"
