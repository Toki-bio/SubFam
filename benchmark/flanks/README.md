# Flank test (2026-11-01, DRAGEN, SubFam fca3b01)

The 3,000 chr19 B2-family copies of make_b2.sh cut again with 0, 25, 100, 500 bp of their own genomic flank on each side (inner part verified identical to the element-only copy for all 3,000), SubFam -n 50. Noise floor: the element-only copies in two shuffled input orders (s1, s2). Compared with the element-only run without labels; the rows also scored against the Dfam and Vassetzky 2021 consensuses (best local identity, coverage >= 0.9). Wall time / memory per run in flanks.tsv. Scripts flanks.sh, flanks2.sh.

## Copies against the rows built from them (copies_vs_rows.sh, 2026-11-01)

Best local identity (coverage of the query >= 0.5) of each of the 3,000 element-only copies and of each of the 60 rows of the flank-0 run to the nearest published consensus. Copies without any hit at that coverage (132 for Dfam, 127 for V2021; short or truncated) are left out; every row has a hit.

| reference | copies: median (IQR), mean | rows: median (IQR), mean |
|---|---|---|
| Dfam B2/B3 (5) | 0.799 (0.751–0.895), 0.818 | 0.915 (0.883–0.964), 0.922 |
| Vassetzky et al. 2021 B2 (5) | 0.817 (0.763–0.895), 0.826 | 0.936 (0.913–0.958), 0.935 |

The consensus of a chunk is closer to the published subfamily consensus than its copies are: the private mutations of single copies drop out.
