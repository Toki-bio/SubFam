# Mouse B2 SINEs (mm39 chr19): SubFam against two published classifications

Run 2026-10-10 in the cloud session (8 threads, MAFFT 7.505, SubFam 1.3.0). Everything is made by
`make_b2.sh WORKDIR` (downloads, extraction, sampling with seed 1, SubFam, grouped alignment).
Result file: `examples/b2_mm39_chr19_n50.aln.fasta` (the five Vassetzky 2021 consensuses B2a–e are its first five rows, `REFERENCE_KRAMEROV_*`; the Dfam ones are likewise in the file grouped by B2a–e).

## Data
- UCSC mm39 `rmsk.txt.gz` and `chr19.fa.gz`; every RepeatMasker copy named B2_Mm1a, B2_Mm1t,
  B2_Mm2, B3 or B3A on chr19, at least 170 bp, without N, on the element strand: 5,972 copies
  (B2_Mm1a 505, B2_Mm1t 541, B2_Mm2 1,384, B3 2,493, B3A 1,049).
- A random 3,000 of them (seed 1; B2_Mm1a 246, B2_Mm1t 274, B2_Mm2 672, B3 1,266, B3A 542), so
  that `-n 50` gives 60 rows.
- Dfam 4.0 consensuses DF000001735–DF000001739 (B2_Mm1a, B2_Mm1t, B2_Mm2, B3, B3A) from the Dfam
  API, as `PUBLISHED_DFAM_*` rows.
- Background on the two classifications: `docs/research/B2_literature.md`. In short, Repbase/Dfam
  names rest on single-nucleotide diagnostics; the indel-based scheme of Vassetzky et al. 2021
  (B2a–e) treats Mm1a + Mm1t + Mm2 as one subfamily, B2a, with Mm1 and Mm2 as "extreme cases in a
  continuum", and B3 as B2d. The B2a–e consensuses exist only in Figure S2 of that paper's
  supplement; `kramerov_consensus.sh` fetches the supplement from Europe PMC and extracts them
  (`kramerov2021_B2_consensus.fa`: B2a 187, B2b 190, B2c 190, B2d 216, B2e 196 bp, IUPAC codes kept).

## Run
`SubFam.sh -n 50 -t 8 -x b2 b2_chr19_3000.fa`: 60 chunk consensuses, 19 s.
`build_examples.py ... --refs dfam_B2.fa --both`: for each rmsk class the consensus of all its
member copies (`CLASSCONS_<class>_consensus_of_<class>_copies`, coverage-relative plurality 0.36), then
the Dfam consensus, then the SubFam rows whose majority is that class, in chunk order.

## What it shows (`rows.tsv`, from `examples/check_rows.py`)
Rows per class (majority): B2_Mm1a 4, B2_Mm1t 6, B2_Mm2 14, B3 27, B3A 9. Identity of the 60 rows
to the head of their own group: median 0.959, min 0.836; 5 rows are closer to another head.

| class | rows | purity of the chunks (p) | identity to own head | reading |
|---|---|---|---|---|
| B2_Mm1a | 4 | 0.78–1.00 | 0.990–1.000 | clean |
| B2_Mm2 | 14 | 0.44–1.00 (11 of 14 ≥ 0.76) | 0.867–1.000; 1 of 14 closer to the Mm1t head | mostly clean |
| B2_Mm1t | 6 | 0.44–0.96 (5 of 6 ≤ 0.66) | 0.957–0.979; 3 of 6 closer to the Mm1a or Mm2 head | no chunk of its own: its copies sit between Mm1a and Mm2 |
| B3 | 27 | 0.56–1.00 | 0.914–0.990 | clean blocks, old (the copies are 23–27 % diverged) |
| B3A | 9 | 0.50–0.92 | 0.836–0.927; 1 of 9 closer to the B3 head | old and mixed with B3 |

The chunk order puts B3/B3A first (chunks 1–21), then the B2 block (22–45: 22–29 Mm2, 30–34 Mm1t/Mm2,
35–38 Mm1a, 39–45 Mm1t/Mm2), then B3 again (46–60): the k-mer tree separates the B2 and B3 lineages;
inside B2 the three Repbase names interleave, with Mm1a as a short pure run and Mm1t never pure.

Read with the literature: this is what the 2021 indel scheme predicts (Mm1a/Mm1t/Mm2 = one
subfamily with a gradient; B3 = another), and not what five separate Repbase entries would
predict for Mm1t. SubFam asserts nothing here; the rows are the data.

## Caveats
- One chromosome, 3,000 of 5,972 copies; the labels are RepeatMasker's (2010 library, see
  `benchmark/line_real/README.md` §2 for the same caveat on hg38).
- B2b, B2c, B2e of the 2021 scheme have no Repbase name and are labelled B2_Mm2/B3/B3A by
  RepeatMasker; the next section labels the copies by the 2021 consensuses instead.
- Purity is against rmsk labels, which are themselves uncertain between Mm1a/Mm1t/Mm2 (2 and 16
  nt apart over 193 bp).

## The same copies against the 2021 subfamilies (B2a–e)

`label_by_consensus.py`: every copy labelled by its best hit to the five 2021 consensuses
(vsearch global alignment, identity threshold 0.5, loosened word prefilter; the way Borodulina et al.
2025 assigned copies). All 3,000 copies get a label (1 tie). Cross-table with the RepeatMasker label:

| rmsk \ 2021 | B2a | B2b | B2c | B2d | B2e |
|---|---|---|---|---|---|
| B2_Mm1a | 246 | 0 | 0 | 0 | 0 |
| B2_Mm1t | 273 | 1 | 0 | 0 | 0 |
| B2_Mm2 | 563 | **104** | 2 | 2 | 1 |
| B3 | 0 | 8 | 141 | **927** | 190 |
| B3A | 2 | 12 | 53 | 248 | 226 |

Median identity of the copies to their consensus: B2a 89.3 %, B2b 81.1 %, B2c 75.1 %, B2d 70.4 %,
B2e 67.2 %. Identity between consensuses (vsearch, global): the three Dfam B2 entries are 94–95 %
to B2a and 90–98 % to each other; Dfam B3 is 92 % to B2d; Dfam B3A is 81 % to B2e and 78 % to B2d
(closest to neither); B2c/B2d/B2e are 77–85 % to each other.

Reading: Mm1a + Mm1t + Mm2 = B2a, as the paper says, except that **104 of 672 rmsk B2_Mm2 copies
are B2b**, a subfamily RepeatMasker has no name for; B3 = B2d plus a share of B2e and B2c; B3A is
not a subfamily in the 2021 scheme but a mixture of old B2d/B2e/B2c copies, and 214 of its 542
copies have no hit at 50 % with the stricter prefilter.

**Grouped by the 2021 label** (`examples/b2_mm39_chr19_n50_kramerov.aln.fasta`, `rows_kramerov.tsv`):
B2a 22 rows (purity 0.50–1.00, 20 of 22 ≥ 0.92), B2b 2 rows (0.48, 0.60), B2c 2 rows (0.92, 0.48),
B2d 34 rows (0.48–0.86), B2e none: its 417 copies never form a chunk majority, they are spread over
the B2d chunks. Identity to own head median 0.951, min 0.829; 12 rows closer to another head: 9 B2d rows
closer to the B2c head, 2 to the B2e head, 1 B2a row to the B2b head.

So the young, indel-defined subfamilies resolve (B2a as a block of 22 pure chunks; B2b and B2c as
their own chunks, small because they are 4–7 % of the copies), and the old ones (B2d, B2e, at 67–70 %
identity to their consensus) do not: their copies are as far from each other as from the other
subfamily's consensus, the best-hit labels are correspondingly uncertain, and the k-mer order mixes
them. That is the limit of any copy-level method on 25–30 % diverged SINEs, and the alignment shows
it rather than hiding it.

**For manual inspection**: `examples/b2_mm39_chr19_n50_both.aln.fasta` has the five Dfam and the
five 2021 consensuses at the top and the 60 chunk consensuses below in chunk order (the k-mer
guide-tree order), each named `<chunk>|rmsk:<label>_p<purity>|kram:<label>_p<purity>`. Things to
look at: the 11-bp central deletion shared by B2a and B2b against B2c/d/e; the 28-bp 3' insertion of
B2d; the τ/terminator region before the A-tail where B2a differs; the run of chunks 22–45 (B2a, from
Mm2-like to Mm1a-like) and where the two B2b chunks (44–45) and the B2c chunks (13–14) sit.

## Baseline: identity-threshold clustering (`../baseline_vsearch/vsearch_sweep.py`, output `vsearch_sweep.txt`, also copied to `../baseline_vsearch/b2_sweep.txt`)

The same 3,000 copies, `vsearch --cluster_fast` (global identity, `--iddef 2`, plus strand) at six thresholds, against
SubFam's 60 chunks and a single family-wide consensus. A group of at least 10 copies counts as a real group; purity
is the fraction of copies carrying the majority label of their group.

| method | groups | singletons | copies in groups ≥ 10 | purity (RepeatMasker label) of those groups |
|---|---|---|---|---|
| one family-wide consensus | 1 | 0 | 100 % | 0.422 |
| SubFam chunks, N = 50 | 60 | 0 | 100 % | 0.793 |
| vsearch id 0.60 | 632 | 389 | 63 % | 0.666 |
| vsearch id 0.70 | 1,613 | 1,488 | 40 % | 0.573 |
| vsearch id 0.80 | 1,992 | 1,906 | 29 % | 0.796 |
| vsearch id 0.90 | 2,560 | 2,474 | 9 % | 0.758 |
| vsearch id 0.95 | 2,868 | 2,812 | 1 % | 0.818 |
| vsearch id 0.98 | 2,976 | 2,960 | 0 % | none |

Reading: these copies are 25–30 % diverged from each other, so a threshold either makes thousands of singletons
(0.80 and above: 1,906 to 2,960 of 3,000 copies alone) or, at 0.60, merges a third of the copies into 43 groups
and leaves the rest alone, with lower purity than SubFam's chunks, which cover every copy. There is no threshold
at which a threshold method gives a short alignment that includes all copies. Caveats: (1) purity against the 2021
labels is not shown in the table because those labels are themselves best hits by the same vsearch identity, so
a vsearch clustering matches them by construction (the full table in `vsearch_sweep.txt` has the column); (2) purity
rewards small groups, which is why the coverage column is printed; (3) `cluster_fast` is greedy and
order-dependent; (4) this is one chromosome and one family, so it illustrates the argument and is not a benchmark.

## Sensitivity (2026-10-11; `sensitivity.sh`, output `sensitivity.tsv`)
Chunk purity of the 3,000-copy run against the RepeatMasker label and the Vassetzky 2021 best-hit label, and the mean `agree` of `chunkstats.tsv`, for other settings and other random samples of the 5,972 copies:

| setting | chunks | purity rmsk | purity 2021 |
|---|---|---|---|
| N = 20 / 30 / 50 / 100 (k 6) | 150 / 100 / 60 / 30 | 0.812 / 0.803 / 0.793 / 0.757 | 0.785 / 0.780 / 0.773 / 0.760 |
| k = 4 / 5 / 6 / 8 / 10 (N 50) | 60 | 0.758 / 0.761 / 0.793 / 0.763 / 0.754 | 0.779 / 0.771 / 0.773 / 0.770 / 0.761 |
| p = 0.25 / 0.36 / 0.5 (N 50, k 6) | 60 | 0.793 (all) | 0.773 (all) |
| MAFFT guide-tree order (`-m`, N 50) | 60 | 0.734 | 0.764 |
| samples seed 2 / 3 / 4 (N 50, k 6) | 60 | 0.778 / 0.781 / 0.774 | 0.781 / 0.777 / 0.779 |

Purity moves by about 0.05 with N (smaller chunks are purer, as expected), by up to 0.04 with k (k = 6 best against rmsk, flat against the 2021 labels), not at all with p (p changes the consensus, not the chunk membership), and by 0.02 between random samples. The MAFFT order is 0.06 lower against rmsk. No setting was tuned on these labels.

## Stress case: 30,000 copies of mm39 chr1–3 (2026-10-11; `stress_chr123/`)
`stress_chr123/pipeline.sh [WORKDIR] [THREADS]` (downloads its own data; chromosomes 1–3 of mm39 from UCSC, the extraction rule of `make_b2.sh`: 41,989 full-length copies,
chr1 14,690, chr2 16,345, chr3 10,954; a random 30,000, seed 1, `subsample.py`), SubFam 1.4.1 for `-n 50` and `-n 50 -P`, 1.4.2 for `-n 20`
(the three runs differ only in the final-alignment mode), 8 threads on a 4-core machine; peak RSS of the whole run. Scored by `score_big.py`
(outputs `score_b2.txt`, `score_b2_n20.txt`): purity against the RepeatMasker class, recovery of the five Dfam consensuses and of the 60 chr19
rows by `../alu_hg38/recover.py` (coverage ≥ 0.9). Consensuses, chunk statistics and chunk tables of the three runs are in `stress_chr123/`;
the 600-row alignment is `examples/b2_mm39_chr123_30k_n50.aln.fasta`.

| run | wall | peak RSS | rows (distinct) | purity | mean agree | Dfam B2_Mm1a / Mm1t / Mm2 / B3 / B3A best identity | chr19 rows found at ≥ 0.95 |
|---|---|---|---|---|---|---|---|
| `-n 50`, k-mer order | 13 min 11 s | 4.1 GB | 600 (588) | 0.805 | 0.795 | 1.000 / 0.984 / 0.995 / 0.943 / 0.889 | 35 of 60 |
| `-n 50 -P` | 7 min 19 s | 0.7 GB | 600 (592) | 0.756 | 0.793 | 1.000 / 0.985 / 0.995 / 0.939 / 0.890 | 38 of 60 |
| `-n 20`, k-mer order (1.4.2) | 15 min 26 s | 4.1 GB | 1,500 (1,486) | 0.826 | 0.801 | 1.000 / 0.985 / 0.995 / 0.943 / 0.890 | 35 of 60 |

Seeds 2 and 3 (`stress_chr123/seeds.sh`, SubFam 1.4.2, `-n 50`, k-mer order, CPU otherwise idle): 6 min 23 s and 6 min 15 s, 4.1 GB; 600 rows (591 and 590 distinct);
purity 0.804 and 0.809; mean agree 0.795 and 0.796; Dfam B2_Mm1a 1.000 / 1.000, Mm1t 0.985 / 0.985, Mm2 0.995 / 0.990, B3 0.954 / 0.952, B3A 0.888 / 0.898; chr19 rows found at ≥ 0.95: 40 of 60 (seed 2).
So the three samples agree to 0.005 in purity and 0.01 in recovery, and the 13 min of seed 1 was the shared-CPU figure (6.3 min on an idle machine).

Reading: at ten times the input the picture of the 3,000-copy run holds (purity 0.805 against 0.793; the three young Dfam consensuses in the rows at
≥ 0.984, the two old ones at 0.89–0.94 as before); `-P` saves most of the memory and half the time at a cost of 0.05 in purity; `-n 20` gives 1,500 rows
at purity 0.826. The old rule for the final alignment (full L-INS-i refinement, single-threaded) took the `-n 20` run past 1.5 h in the refinement
alone before it was stopped; 1.4.2 aligns the 1,500 rows with MAFFT `--auto` in seconds. The 11.5 min of the `-n 20` ordering step was measured
while another process shared the CPU. The chr19 rows found among the 600 at ≥ 0.95 are mostly the B2 rows; the B3/B3A rows of the two runs differ,
as expected for 25–30 % diverged copies drawn from different chromosomes.
