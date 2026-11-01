# Flank test (2026-11-01, DRAGEN, SubFam fca3b01)

The 3,000 chr19 B2-family copies of make_b2.sh cut again with 0, 25, 100, 500 bp of their own genomic flank on each side (inner part verified identical to the element-only copy for all 3,000), SubFam -n 50. Noise floor: the element-only copies in two shuffled input orders (s1, s2). Compared with the element-only run without labels; the rows also scored against the Dfam and Vassetzky 2021 consensuses (best local identity, coverage >= 0.9). Wall time / memory per run in flanks.tsv. Scripts flanks.sh, flanks2.sh.

## Copies against the rows built from them (copies_vs_rows.sh, 2026-11-01)

Best local identity (coverage of the query >= 0.5) of each of the 3,000 element-only copies and of each of the 60 rows of the flank-0 run to the nearest published consensus. Copies without any hit at that coverage (132 for Dfam, 127 for V2021; short or truncated) are left out; every row has a hit.

| reference | copies: median (IQR), mean | rows: median (IQR), mean |
|---|---|---|
| Dfam B2/B3 (5) | 0.799 (0.751–0.895), 0.818 | 0.915 (0.883–0.964), 0.922 |
| Vassetzky et al. 2021 B2 (5) | 0.817 (0.763–0.895), 0.826 | 0.936 (0.913–0.958), 0.935 |

The consensus of a chunk is closer to the published subfamily consensus than its copies are: the private mutations of single copies drop out.

## Ten input permutations and N runs (permutations_and_N.sh, Monsoon, 2026-11-01)

Ten shuffles (seeds 101-110) of the 3,000 clean B2 copies: original rows recurring at >= 0.95: 37-40/60 (median 38); median best identity of an original row 0.96; copy pairs kept in the same chunk 0.35-0.58 (median 0.39). With the two earlier shuffles (47, 36) the range over twelve orders is 36-47/60.
N runs: 1,000 of the 3,000 copies carry one interior N run of 1-10 bases. SubFam 1.5.1 (k-mers span the N) and the same code with the k-mer window restarted at an invalid letter: 37/60 original rows recur in both, pairs kept against the N-free run 0.166 and 0.170, against each other 0.195. N runs (which also replace 1-10 bases of 1,000 copies) perturb chunk membership more than a reshuffle does (0.17 against a median of 0.39); restarting the k-mer window at an N does not reduce this, and is not worth a release before the paper.
