# Flank test (2026-11-01, DRAGEN, SubFam fca3b01)

The 3,000 chr19 B2-family copies of make_b2.sh cut again with 0, 25, 100, 500 bp of their own genomic flank on each side (inner part verified identical to the element-only copy for all 3,000), SubFam -n 50. Noise floor: the element-only copies in two shuffled input orders (s1, s2). Compared with the element-only run without labels; the rows also scored against the Dfam and Vassetzky 2021 consensuses (best local identity, coverage >= 0.9). Wall time / memory per run in flanks.tsv. Scripts flanks.sh, flanks2.sh.
