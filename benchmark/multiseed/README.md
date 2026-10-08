# Eight simulated replicates per scenario (2026-10-08)

`run.sh` repeats the README benchmark for seeds 1-8 (`simulate.py --seed`), SubFam `-n 50` / `-n 20` against VSEARCH `--cluster_fast --consout` at seven identity thresholds,
scored with `evaluate.py` ("recovered" = subfamily source matched at >= 99 % identity by a row whose copies are mostly that subfamily). `results.tsv` has every run.
Seed 1 reproduces the numbers in the main README (checked for the middle scenario: SubFam -n 50 40 rows 8/8, -n 20 100 rows 8/8, VSEARCH 0.80 210 rows 8/8, 0.85 1,133 rows 8/8).

| scenario | SubFam -n 20 | SubFam -n 50 | VSEARCH: best single --id over the 8 seeds | VSEARCH: oracle id chosen per seed |
|---|---|---|---|---|
| young (3 %) | 8.00 recovered in 8/8 runs, 100 rows | 8.00, 40 rows | 0.95: 8.00 (1,278 rows); 0.90: 8.00 (92 rows) | 8.00, median 92 rows |
| middle (8 %) | 8.00 in 8/8 runs, 100 rows | 7.88 (7-8), 40 rows | 0.85: 8.00 (1,119 rows); 0.80: 7.88 (194 rows) | 8.00, median 194 rows |
| old (15 %) | 8.00 in 8/8 runs, 100 rows | 7.38 (6-8), 40 rows | 0.70: 5.25 (3-7) with 538 rows | 5.25, median 529 rows |

Reading: the identity that works for VSEARCH moves with the age of the family (0.90-0.95, 0.80-0.85, 0.70); SubFam -n 20 recovered all eight sources in all 24 runs with 100 rows.
`-n 50` loses the smallest subfamilies (40 and 60 copies, below N) in some old and middle runs: N is the resolution parameter, in copies.
Limits: the simulation is star-like (each copy decays independently from its subfamily master), balanced tree, 2.7 % between sister subfamilies, no CpG effect, no recombination,
no truncation, one family length (300 bp). "Oracle id" = for each seed the threshold with the most subfamilies recovered, ties broken by fewest rows (uses the truth, so it is the best case for VSEARCH).
