# Baseline: MeShClust2 (mean-shift, alignment-free identity) on the three real sets (2026-10-11)

MeShClust2 (Girgis lab, github.com/BioinformaticsToolsmith/MeShClust2, commit d7fdb42; built with g++ 13.3 after adding
`-include cstdint -include limits -include algorithm` to its `CMAKE_CXX_FLAGS`), run at several `--id` values with 8 threads on the
same copies as `../baseline_vsearch`, scored by the same rule (`meshclust_sweep.py`: purity = share of copies carrying the majority
label of their group; a real group has ≥ 10 copies). Raw outputs: `alu_mesh.txt`, `b2_mesh.txt`, `l1_mesh.txt`; driver `run_meshclust.sh`.
MeShClust2 was chosen as the threshold-free-in-spirit comparator (its `--id` is a cutoff, but assignment is by mean shift, not greedy).
MeShClust v3 (the Identity repository) was built too and not run.

| set | method | groups | singletons | copies in groups ≥ 10 | purity of those | wall s |
|---|---|---|---|---|---|---|
| B2, 3,000 | **SubFam N = 50** | **60** | **0** | **100 %** | **0.79** | 20 |
| | meshclust2 0.60 | 636 | 533 | 73 % | 0.46 | 9 |
| | meshclust2 0.70 | 2,021 | 1,944 | 29 % | 0.53 | 10 |
| | meshclust2 0.80 | 2,328 | 2,233 | 16 % | 0.56 | 9 |
| | meshclust2 0.90 / 0.95 | 2,738 / 2,882 | 2,711 / 2,865 | 7 / 3 % | 0.74 / 0.95 | 12 / 11 |
| L1PA2–8, 2,222 | **SubFam N = 20, -c** | **111** | **0** | **100 %** | **0.70** | 391 |
| | meshclust2 0.80 | 317 | 246 | 81 % | 0.27 | 8 |
| | meshclust2 0.90 | 709 | 639 | 66 % | 0.34 | 9 |
| | meshclust2 0.95 | 1,992 | 1,957 | 7 % | 0.65 | 11 |
| | meshclust2 0.98 | 2,219 | 2,216 | 0 % | – | 11 |
| Alu loci, 316 | **SubFam N = 20** | **16** | **0** | **100 %** | **0.92** | 2 |
| | meshclust2 0.90 | 1 | 0 | 100 % | 0.52 | 4 |
| | meshclust2 0.95 | 11 | 4 | 95 % | 0.77 | 5 |
| | meshclust2 0.98 | 40 | 24 | 78 % | 0.86 | 4 |
| | meshclust2 0.99 | 73 | 61 | 72 % | 0.96 | 9 |

Reading: the same shape as vsearch. On the old families no `--id` gives a short, complete result: at 0.60–0.80 the real groups are impure
(B2 0.46–0.56, L1 0.27–0.34, against SubFam's 0.79 and 0.70), at 0.90 and above most copies are singletons. On the young Alu loci
MeShClust2 at 0.95 gives 11 groups covering 95 % of the copies at purity 0.77 (vsearch 0.95: 18 groups, 89 %, 0.94; SubFam: 16 rows, 100 %, 0.92),
and at 0.90 one cluster. MeShClust2 is much faster than SubFam on L1 (8–11 s against 6.5 min) because it does not align or build consensuses.
Caveats as for vsearch: labels are annotations, one chromosome per family, purity rewards small groups. Its 2018 and v3 versions were not run.
