# Baseline: identity-threshold clustering against SubFam's chunks (three real-data examples)

`vsearch --cluster_fast` (vsearch 2.27.0; global identity, `--iddef 2`, plus strand) at several identities, scored
against SubFam's chunks and against a single family-wide consensus. A "real" group has at least 10 copies. Purity =
fraction of copies carrying the majority label of their group. Scripts: `vsearch_sweep.py`, `run_sweeps.sh`
(L1, Alu); B2 is step 6 of `../sine_b2/make_b2.sh`. Outputs: `b2_sweep.txt`, `l1_sweep.txt`, `alu_sweep.txt`.
Run 2026-10-10, SubFam 1.3.0, 8 threads. L1 SubFam run: 6 min 31 s (its own, `-n 20 -c`).

| example | copies | divergence between copies | method | groups | singletons | copies in groups ≥ 10 | purity of those groups | purity of all groups |
|---|---|---|---|---|---|---|---|---|
| B2 (rmsk label) | 3,000 | 25–30 % | one consensus | 1 | 0 | 100 % | 0.42 | 0.42 |
| | | | **SubFam, N = 50** | **60** | **0** | **100 %** | **0.79** | **0.79** |
| | | | vsearch 0.60 | 632 | 389 | 63 % | 0.67 | 0.76 |
| | | | vsearch 0.70 | 1,613 | 1,488 | 40 % | 0.57 | 0.83 |
| | | | vsearch 0.80 | 1,992 | 1,906 | 29 % | 0.80 | 0.94 |
| | | | vsearch 0.90 | 2,560 | 2,474 | 9 % | 0.76 | 0.97 |
| | | | vsearch 0.95 | 2,868 | 2,812 | 1 % | 0.82 | 0.995 |
| L1PA2–8 (rmsk label) | 2,222 | 1–7 % from own consensus; ≥ 1 kb, 5′-truncated | one consensus | 1 | 0 | 100 % | 0.22 | 0.22 |
| | | | **SubFam, N = 20, `-c`** | **111** | **0** | **100 %** | **0.70** | **0.70** |
| | | | vsearch 0.80 | 45 | 31 | 96 % | 0.36 | 0.38 |
| | | | vsearch 0.90 | 753 | 686 | 60 % | 0.48 | 0.68 |
| | | | vsearch 0.95 | 1,557 | 1,475 | 22 % | 0.83 | 0.95 |
| | | | vsearch 0.98 / 0.99 | 2,178 / 2,217 | 2,152 / 2,212 | 0 % | none | 1.00 |
| Alu, Konkel 2015 loci (Price label) | 316 | young (AluY, Ya5, Yb8) | one consensus | 1 | 0 | 100 % | 0.52 | 0.52 |
| | | | **SubFam, N = 20** | **16** | **0** | **100 %** | **0.92** | **0.92** |
| | | | vsearch 0.95 | 18 | 8 | 89 % | 0.94 | 0.95 |
| | | | vsearch 0.97 | 34 | 15 | 78 % | 0.96 | 0.97 |
| | | | vsearch 0.98 | 70 | 38 | 65 % | 0.97 | 0.98 |
| | | | vsearch 0.99 | 140 | 101 | 35 % | 0.98 | 0.99 |

## Reading
- **Old, diverged families (B2 at 25–30 %, L1 over 40 Myr).** There is no threshold with a short result that keeps all copies. Low thresholds
  merge across classes (L1 at 0.80: 45 groups, three of them hold 96 % of the copies, purity 0.36, below SubFam's 0.70), high thresholds
  leave most copies alone (B2 at 0.80: 1,906 of 3,000 copies are singletons). At no setting is the output both short and complete.
- **Young family (Konkel Alu).** The threshold method does as well as SubFam: vsearch at 0.95 gives 18 groups covering 89 % of the copies at purity 0.94,
  against SubFam's 16 chunks at 0.92 covering all. For young sequences where a threshold exists, SubFam's gain is completeness
  (no singletons) and the ordered alignment, not purity.
- This comparison favours neither method in general: purity rewards small groups (the coverage column corrects for that); SubFam's chunk
  size is a free parameter, set before looking at labels (N = 50 / 20 / 20, the repository defaults for the three sizes); labels are RepeatMasker's
  or the Price et al. subfamily, not truth; `cluster_fast` is greedy and order-dependent; one chromosome per family.
- Not run: CD-HIT (same logic as vsearch), COSEG (needs the 5′-complete copies and a different input), MeShClust (not installed here).
