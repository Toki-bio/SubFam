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

## How faithfully do rows stand for their copies? (seeds 1-4, 8,000 copies per scenario)
`fidelity.py`: every copy is aligned (vsearch --usearch_global, --iddef 1) against all SubFam consensuses of its run. "Own row is best" = the row built from the copy's own chunk gives its highest identity.
Ceiling = the same test against the 8 true subfamily masters (is the copy nearest its true master).

| scenario | -n | own row is the copy's best row | median identity copy -> own row | 5th percentile | ceiling (true master nearest) |
|---|---|---|---|---|---|
| young (3 %) | 50 / 20 | 98.1 % / 95.0 % | 96.7 % | 94.1 % | 100.0 % |
| middle (8 %) | 50 / 20 | 87.3 % / 88.0 % | 91.0 % | 86.2 % | 100.0 % |
| old (15 %) | 50 / 20 | 57.5 % / 48.9 % | 82.6 % | 76.3 % | 99.8 % |

Reading: the rows are faithful summaries for young and middle families and much less so for old ones, where chunks mix sister subfamilies (chunk purity 0.71 in the README), even though every row still
contains its source consensus. A copy that is not nearest to its own row is still represented in `chunks.tsv`; the claim is organisation, not per-copy attribution.

Re-run 2026-10-10 with SubFam 1.3.0 (`run.sh`, 8 threads): all 216 rows of `results.tsv` are byte-identical to the file above.
