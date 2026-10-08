# Harder simulation (CpG decay, source elements, age mixture, truncation)

`simulate2.py` extends `benchmark/simulate.py`: a CpG-rich ancestor with CpG decay at 6 x the rate of other sites (C>T, G>A); 1-5 active sources per subfamily with skewed activity (the dominant source is the subfamily master, the others carry private mutations);
per-copy age drawn around the scenario mean (gamma, CV 0.35); transitions twice as likely as transversions; 10 % indels; a fraction of copies 5'-truncated (100 bp + exponential, mean 90 bp).
`evaluate2.py` scores representatives against the true subfamily masters from a local alignment, with coverage (>= 0.95 of the master) and two identities: over all columns, and leaving out the master's CpG positions (decayed copies make every consensus TG/CA there, so the consensus is not expected to equal the master at CpG).
`run_hard.sh` runs 5 scenarios x seeds; `results.tsv` = seeds 1-4. All parameters are mine, chosen to be harder than the simple model, not calibrated to a genome.

Scenarios (8 subfamilies, 2,000 copies, 300 bp): young / middle / old = mean divergence 3 / 8 / 15 % with 30 % of copies truncated; mixed = four young and four old subfamilies in one dataset; heavytrunc = middle with 80 % truncated copies.

Subfamilies recovered of 8 (mean over 4 seeds; "noCpG" = identity >= 99 % over non-CpG columns with coverage >= 0.95 and a group mostly that subfamily; "all" = same over all columns), rows = median:

| scenario | method | rows | recovered, all columns | recovered, noCpG |
|---|---|---|---|---|
| young | SubFam -n 50 / -n 20 | 40 / 100 | 7.25 / 7.75 | 7.25 / 7.75 |
| | VSEARCH 0.95 / 0.90 / 0.85 / 0.80 | 1,808 / 965 / 283 / 50 | 4.5 / 7.25 / 7.5 / 7.0 | 5.0 / 7.75 / 7.75 / 7.25 |
| middle | SubFam -n 50 / -n 20 | 40 / 100 | 4.5 / 4.75 | 6.5 / 7.5 |
| | VSEARCH 0.85 / 0.80 / 0.75 / 0.70 | 1,696 / 1,067 / 369 / 75 | 1.25 / 0.75 / 0.5 / 0.25 | 6.25 / 7.5 / 7.0 / 6.25 |
| old | SubFam -n 50 / -n 20 | 40 / 100 | 0 / 0 | 4.75 / 6.25 |
| | VSEARCH 0.80 / 0.75 / 0.70 | 1,724 / 999 / 385 | 0 / 0 / 0 | 3.5 / 5.5 / 4.5 |
| mixed ages | SubFam -n 50 / -n 20 | 40 / 100 | 4.0 / 4.0 | 6.25 / 7.0 |
| | VSEARCH 0.90 / 0.85 / 0.80 / 0.75 | 1,404 / 973 / 769 / 451 | 3.75 / 4.0 / 4.0 / 3.25 | 4.0 / 4.0 / 4.25 / 6.0 |
| 80 % truncated | SubFam -n 20 / -n 20 -c / -m -n 20 | 100 | 1.5 / 1.5 / 0.0 | 6.5 / **7.25** / 0.75 |
| | VSEARCH 0.85 / 0.80 / 0.75 | 1,597 / 871 / 295 | 0 / 0 / 0 | 3.0 / 4.75 / 3.0 |

What it shows:
- Not 8/8 any more: with CpG decay, source elements, age mixture and truncation SubFam -n 20 recovers 6.25-7.75 of 8 (noCpG), still at least as many as VSEARCH at its best fixed identity in every scenario (equal in young and middle, 6.25 vs 5.5 old, 7.0 vs 6.0 mixed, 7.25 vs 4.75 heavy truncation with -c) and with 100 rows against 50-1,700.
- CpG: identity to the master over all columns is below 99 % for essentially every consensus (old: 0 recovered by any method), because CpG decay moves the consensus; the noCpG column is the fair one, both are reported.
- Truncation: with 30 % truncated copies -c changes nothing (the default plurality still covers the full length); with 80 % it lifts recovery 6.5 -> 7.25, and MAFFT guide-tree ordering (-m) collapses (0.75), so the k-mer tree is the better ordering for truncated copies.
- Purity of the chunks falls to 0.53-0.87 (mean over seeds), see results.tsv: the ordering is harder to get right than in the simple model.
Limits: 4 seeds, parameters not calibrated to real data; one length (300 bp); "oracle" identity thresholds not tested here (the VSEARCH rows are fixed identities, best one shown per scenario in the paper table should be chosen the same way as in benchmark/multiseed).
