#!/usr/bin/env bash
# The threshold baseline on the three real-data examples: vsearch --cluster_fast (global identity, --iddef 2,
# plus strand) at several identities, against SubFam's chunks and one family-wide consensus.
# B2 is produced by benchmark/sine_b2/make_b2.sh (step 6). This script does L1 and Alu from the inputs named in
# examples/rebuild.sh ($L1 = l1.fa of benchmark/line_real/extract_l1.py, $KONKEL = alu_bodies.fa of benchmark/alu_konkel/prep.py).
#   run_sweeps.sh WORKDIR L1.fa ALU_BODIES.fa [THREADS]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$HERE/../..
W=$1; L1=$(readlink -f "$2"); ALU=$(readlink -f "$3"); T=${4:-4}; mkdir -p "$W"; cd "$W"
# L1PA2-8, chr4 hg38: SubFam -n 20 -c (as examples/rebuild.sh), label = RepeatMasker class (the first field of the id)
grep '^>' "$L1" | sed 's/^>//' | awk -F_ '{ print $0 "\t" $1 }' > l1_truth.tsv
"$TOP/SubFam.sh" -n 20 -c -t "$T" -o sf_l1 -x l1 "$L1"
python3 "$HERE/vsearch_sweep.py" "$L1" sf_l1/l1.chunks.tsv l1_truth.tsv rmsk --ids 0.80,0.90,0.95,0.98,0.99 --threads "$T" | tee l1_sweep.txt
# Alu, Konkel et al. 2015 loci: SubFam -n 20, label = the subfamily in the id (acc|label)
awk -F'|' '/^>/ { print substr($1, 2) "|" $2 "\t" $2 }' "$ALU" > alu_truth.tsv
"$TOP/SubFam.sh" -n 20 -t "$T" -o sf_alu -x alu "$ALU"
python3 "$HERE/vsearch_sweep.py" "$ALU" sf_alu/alu.chunks.tsv alu_truth.tsv Price --ids 0.95,0.97,0.98,0.99,0.995 --threads "$T" | tee alu_sweep.txt
