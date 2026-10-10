#!/usr/bin/env bash
# Mouse B2 example: full-length B2-family copies of mm39 chr19, SubFam -n 50, grouped by the
# RepeatMasker label, with the Dfam consensus beside each group head.
#   make_b2.sh WORKDIR [THREADS]
# Needs: curl, python3 (numpy), mafft; SubFam.sh, benchmark/line_real/extract_l1.py and
# examples/build_examples.py of this repository. Everything is deterministic (sample seed 1).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$HERE/../..
W=$1; T=${2:-4}; mkdir -p "$W/dl"; cd "$W"

# 1. data: UCSC mm39 RepeatMasker table and chr19; Dfam 4.0 consensuses of the five B2-family entries
[ -s dl/rmsk.txt.gz ]  || curl -sSf -o dl/rmsk.txt.gz  https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz
[ -s dl/chr19.fa.gz ] || curl -sSf -o dl/chr19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr19.fa.gz
if [ ! -s dfam_B2.fa ]; then
    for a in DF000001735 DF000001736 DF000001737 DF000001738 DF000001739; do   # B2_Mm1a B2_Mm1t B2_Mm2 B3 B3A
        curl -sSf "https://www.dfam.org/api/families/$a/sequence?format=fasta"; echo
    done | grep -v '^$' | sed 's/^>\(DF[0-9.]*\) \(.*\)/>DFAM_\2/' > dfam_B2.fa
fi

# 2. copies: every rmsk B2_Mm1a / B2_Mm1t / B2_Mm2 / B3 / B3A copy on chr19 of at least 170 bp,
#    element strand, no N (5,972 copies); then a random 3,000 (seed 1) so that -n 50 gives 60 rows
python3 -I "$TOP/benchmark/line_real/extract_l1.py" dl/rmsk.txt.gz dl/chr19.fa.gz b2_chr19.fa 170 '^(B2_Mm1a|B2_Mm1t|B2_Mm2|B3|B3A)$'
python3 -I "$HERE/subsample.py" b2_chr19.fa 3000 1 b2_chr19_3000.fa b2_truth.tsv

# 2b. the five mouse B2 subfamily consensuses of Vassetzky et al. 2021 (Figure S2 of the supplement)
bash "$HERE/kramerov_consensus.sh" kramerov2021_B2_consensus.fa supp

# 3. SubFam, then the grouped alignment (member consensus of each rmsk class + Dfam consensus as heads)
"$TOP/SubFam.sh" -n 50 -t "$T" -o sf50 -x b2 b2_chr19_3000.fa | tee sf50.log
python3 "$TOP/examples/build_examples.py" b2_mm39_chr19_n50.aln.fasta b2_chr19_3000.fa b2_truth.tsv \
    sf50/b2.cons.fasta sf50/b2.chunks.tsv --refs dfam_B2.fa --both --extra-heads kramerov2021_B2_consensus.fa --threads "$T" --tmp build_tmp
python3 "$TOP/examples/check_rows.py" b2_mm39_chr19_n50.aln.fasta > rows.tsv

# 4. every copy labelled by its best hit to those consensuses, the same alignment grouped by that label, and one
#    alignment with both labels on every row for manual inspection
python3 "$HERE/label_by_consensus.py" b2_chr19_3000.fa kramerov2021_B2_consensus.fa b2_truth_kramerov.tsv b2_truth.tsv --threads "$T" 2>&1 | tee label.log
python3 "$TOP/examples/build_examples.py" b2_mm39_chr19_n50_kramerov.aln.fasta b2_chr19_3000.fa b2_truth_kramerov.tsv \
    sf50/b2.cons.fasta sf50/b2.chunks.tsv --refs kramerov2021_B2_consensus.fa --both --extra-heads dfam_B2.fa --threads "$T" --tmp build_tmp_k
python3 "$TOP/examples/check_rows.py" b2_mm39_chr19_n50_kramerov.aln.fasta > rows_kramerov.tsv
python3 "$HERE/combine.py" b2_mm39_chr19_n50_both.aln.fasta b2_chr19_3000.fa sf50/b2.chunks.tsv sf50/b2.cons.fasta \
    b2_truth.tsv rmsk b2_truth_kramerov.tsv kram dfam_B2.fa kramerov2021_B2_consensus.fa --threads "$T" --tmp build_tmp_c

# 5. the chaos-to-order figure: 390 hand-picked real copies (B2a 120, B2b 30, B2c 120, B2d 120), SubFam -n 20,
#    random order / k-mer order / chunk consensuses on one alignment (needs matplotlib)
python3 "$HERE/chaos_subset.py" b2_chr19_3000.fa b2_truth_kramerov.tsv chaos_subset.fa chaos_subset.tsv
"$TOP/SubFam.sh" -n 20 -t "$T" -o sf_chaos -x sub chaos_subset.fa
python3 "$HERE/figure_chaos_to_order.py" chaos_subset.fa chaos_subset.tsv sf_chaos/sub.chunks.tsv sf_chaos/sub.cons.fasta chaos_to_order_b2.png --threads "$T"

# 6. the threshold baseline: vsearch --cluster_fast at 0.60-0.98 against SubFam's 60 chunks and one family consensus
python3 "$TOP/benchmark/baseline_vsearch/vsearch_sweep.py" b2_chr19_3000.fa sf50/b2.chunks.tsv b2_truth.tsv rmsk b2_truth_kramerov.tsv V2021 --threads "$T" | tee vsearch_sweep.txt
