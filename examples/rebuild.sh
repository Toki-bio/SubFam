#!/usr/bin/env bash
# Rebuild every alignment in examples/ from scratch, with the exact commands and seeds used.
#   examples/rebuild.sh WORKDIR [THREADS]
# Deterministic given the same MAFFT version (7.505 was used) and SubFam version (1.3.0).
#
# Inputs that are not generated here (see the benchmark READMEs for how to make them):
#   $KONKEL  alu_bodies.fa from benchmark/alu_konkel/prep.py (GenBank KT305395-KT305737 + the
#            Price et al. 2004 supplement; neither is stored in the repository)
#   $L1      l1.fa from benchmark/line_real/extract_l1.py rmsk.txt.gz chr4.fa.gz l1.fa 1000
#            (UCSC hg38; 2,222 L1PA2-8 copies >= 1 kb on chr4)
# The B2 example is built by benchmark/sine_b2/make_b2.sh (downloads its own data).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$HERE/..; SUB=$TOP/SubFam.sh; BUILD=$HERE/build_examples.py
W=$1; T=${2:-4}; KONKEL=${KONKEL:-}; L1=${L1:-}
mkdir -p "$W"; cd "$W"

# 1. simple simulation, middle family (8 %), seed 1, SubFam -n 50
python3 "$TOP/benchmark/simulate.py" middle --divergence 0.08 --seed 1
"$SUB" -n 50 -t "$T" -o sf_simple -x simple middle.fasta
python3 "$BUILD" sim_simple_middle_n50.aln.fasta middle.fasta middle.truth.tsv \
    sf_simple/simple.cons.fasta sf_simple/simple.chunks.tsv --threads "$T" --tmp tmp_simple

# 2. harder simulation (CpG decay, sources, truncation), old family (15 %), seed 1, SubFam -n 20
python3 "$TOP/benchmark/hard/simulate2.py" old --divergence 0.15 --seed 1
"$SUB" -n 20 -t "$T" -o sf_old -x old old.fasta
python3 "$BUILD" sim_hard_old_n20.aln.fasta old.fasta old.truth.tsv \
    sf_old/old.cons.fasta sf_old/old.chunks.tsv --threads "$T" --tmp tmp_old

# 3. Konkel et al. 2015 Alu bodies, SubFam -n 20, Price et al. 2004 consensuses as heads
if [ -s "$KONKEL" ]; then
    awk -F'|' '/^>/ { print substr($1, 2) "|" $2 "\t" $2 }' "$KONKEL" > konkel_truth.tsv   # id = acc|label
    "$SUB" -n 20 -t "$T" -o sf_konkel -x konkel "$KONKEL"
    python3 "$BUILD" konkel_alu_n20.aln.fasta "$KONKEL" konkel_truth.tsv \
        sf_konkel/konkel.cons.fasta sf_konkel/konkel.chunks.tsv \
        --refs "$HERE/refs/price2004_alu_consensus.fa" --threads "$T" --tmp tmp_konkel
fi

# 4. real L1PA2-8 copies of hg38 chr4, SubFam -n 20 -c, grouped by the RepeatMasker label;
#    each head is the consensus of at most 60 member copies (random, seed 1)
if [ -s "$L1" ]; then
    grep '^>' "$L1" | sed 's/^>//' | awk -F_ '{ print $0 "\t" $1 }' > l1_truth.tsv
    "$SUB" -n 20 -c -t "$T" -o sf_l1 -x l1 "$L1"
    python3 "$BUILD" real_l1pa_chr4_n20c.aln.fasta "$L1" l1_truth.tsv \
        sf_l1/l1.cons.fasta sf_l1/l1.chunks.tsv --sample 60 --threads "$T" --tmp tmp_l1
fi

for f in *.aln.fasta; do python3 "$HERE/check_rows.py" "$f" > "${f%.aln.fasta}.rows.tsv"; done
