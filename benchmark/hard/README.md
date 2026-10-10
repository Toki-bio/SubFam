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

## L1-like simulation (`simulate_line.py`, `evaluate_line.py`): first pass, seed 1 (2026-10-08)
6 kb master (CpG-rich 5' 900 bp), 8 subfamilies (20 diagnostic substitutions per branch, half in the 5' region, half in the last 1.5 kb), copies 10 % diverged with CpG decay and source elements,
97 % of the copies 5'-truncated (200 bp + exponential, mean 900), 15 % with a 5'-inverted part, 10 % with an internal deletion of 50-500 bp. Seed 1: 2,000 copies, median length 813 bp, 69 full length, 304 copies of 2 kb or more.
`evaluate_line.py` aligns each row to the masters (vsearch, terminal gaps ignored), identity over the aligned span leaving out master CpG positions, coverage = master positions spanned / master length;
recovered at coverage C = a row from >= 3 copies, mostly one subfamily, identity >= 99 %, covering >= C of the master. (An earlier version of the script scored coverage 1.0 and dropped most rows because vsearch --iddef 1 counts terminal gaps; both fixed before the numbers below.)

| method | rows | chunk purity | median row coverage | recovered at coverage 0.2 / 0.5 / 0.9 |
|---|---|---|---|---|
| SubFam -n 50 | 40 | 0.354 | 0.13 | 1 / 1 / 1 |
| SubFam -n 20 | 100 | 0.451 | 0.14 | 2 / 1 / 1 |
| SubFam -n 20 -c | 100 | 0.451 | 0.16 | **5** / 3 / 2 |
| SubFam -n 20 -k 8 | 100 | 0.461 | 0.14 | 2 / 2 / 2 |
| SubFam -n 20 -m (MAFFT order) | 100 | 0.452 | 0.15 | 0 / 0 / 0 |
| recipe: long copies (>= 2 kb, 304 copies) only, -n 10 -c | 30 | 0.408 (long copies) | 0.52 | 3 / 2 / 2 |
| recipe: long copies only, -n 20 -c | 15 | 0.336 | 0.54 | 4 / 2 / 1 |
| VSEARCH 0.90 / 0.85 / 0.80 / 0.75 (--consout) | 1,984 / 1,790 / 1,221 / 569 | 0.999 / 0.973 / 0.875 / 0.652 | 0.00-0.13 | 0-1 / 0 / 0 |

Placing every copy by best hit (iddef 2) in the rows of the long-only run (n=10) puts 41 % of copies under 500 bp, 54 % of 500-1,000 bp and 56 % of 1-2 kb in a row whose majority is their own subfamily (8 subfamilies: chance is about 12-25 %).
Reading (one seed, so exploratory): truncated L1-like copies are hard for every method. Chunk purity stays at 0.35-0.46 because length confounds the k-mer order (SCALING.md section 4 saw the same). `-c` is what makes the rows reach beyond the 3' end (5 versus 2 recovered at 0.2 coverage)
and the k-mer order is far better than the MAFFT order here too. Rows cover 13-16 % of the master by default and about half of it from the long copies. Identity-threshold clustering recovers essentially none, with 570-1,980 rows.

Seeds 1-3 (new simulations of the same design; copies >= 2 kb: 304, 329, 327; median length 813, 838, 839 bp). Recovered subfamilies of 8 at coverage 0.2 / 0.5 / 0.9, per seed and mean:

| method | seed 1 | seed 2 | seed 3 | mean |
|---|---|---|---|---|
| SubFam -n 20 | 2 / 1 / 1 | 3 / 2 / 2 | 2 / 1 / 1 | 2.3 / 1.3 / 1.3 |
| SubFam -n 20 -c | 5 / 3 / 2 | 5 / 2 / 2 | 4 / 1 / 1 | **4.7** / 2.0 / 1.7 |
| long copies only (>= 2 kb), -n 20 -c (15-16 rows) | 4 / 2 / 1 | 4 / 2 / 2 | 3 / 3 / 2 | 3.7 / 2.3 / 1.7 |

Chunk purity of the 100-row runs: 0.451, 0.467, 0.432. Seed 1 VSEARCH (0.90 to 0.75, 570-1,984 rows) recovered 0-1 of 8 at every coverage; VSEARCH was not run on seeds 2 and 3 (about ten minutes per setting on 6 kb copies, and the batch job was lost before it got there).
So the seed 1 picture holds in three simulations: `-c` roughly doubles the number of sources whose 3' part is recovered (4.7 against 2.3 at coverage 0.2), but the rows cover only about 15 % of the master (half from the long copies only), chunk purity stays at 0.43-0.47, and few sources are recovered over their full length (1.3-2.3 of 8 at 0.9).

## Per-chunk quality statistics (`chunk_stats.py`): prototype, one simulation
Question: can a chunk's name or table carry what kind of chunk it is (one lineage, two lineages, one lineage plus junk) and how good its alignment is?
`chunk_stats.py` reads the chunk alignments kept by `-K` and reports per chunk: n; `id_med` (median identity of members to the chunk's plurality sequence); `agree` (mean fraction of members carrying the plurality base in core columns);
`ragged` (core columns where >= 20 % of members have a gap); a 2-means split of the members on variable columns (`f2` = smaller part, `best` = strongest separating column); `outl` (fraction of members far below the median identity); and a type T1 / T2 / TJ.
Test: harder-simulation middle family, seed 1, `-n 20`, 100 chunks, true composition from the simulation (pure: top subfamily >= 90 %; two: top two each >= 20 % and together >= 90 %; other):

| true composition | chunks | median agree | median id_med | called T1 / T2 / TJ |
|---|---|---|---|---|
| pure | 26 | 0.882 | 0.884 | 16 / 8 / 2 |
| two subfamilies | 43 | 0.851 | 0.852 | 9 / 32 / 2 |
| other (three or more, or diffuse) | 31 | 0.834 | 0.835 | 6 / 23 / 2 |

Separating pure from mixed chunks: AUC 0.92 by `agree`, 0.90 by `id_med`, 0.76 by `best`, 0.60 by `outl`, 0.39 by `f2`. So a simple cohesion number carries most of the signal; the crude split and junk heuristics add little, and the T1/T2/TJ call as written is not reliable
(38 % of pure chunks are called as split). Caveats: one simulation and one seed; cohesion depends on family age, so a threshold has to be relative to the run (for example a chunk's `agree` against the run's median), not absolute; thresholds untuned; no real junk was simulated (junk here = diffuse or old copies).

Re-run 2026-10-10 with SubFam 1.3.0 (`run_hard.sh OUT "1 2 3 4" 8`): all 200 rows of `results.tsv` are byte-identical to the committed file.
