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
member copies (`TRUE_<class>_consensus_of_<class>_copies`, coverage-relative plurality 0.36), then
the Dfam consensus, then the SubFam rows whose majority is that class, in chunk order.

## What it shows (`rows.tsv`, from `examples/check_rows.py`)
Rows per class (majority): B2_Mm1a 4, B2_Mm1t 6, B2_Mm2 14, B3 27, B3A 9. Identity of the 60 rows
to the head of their own group: median 0.959, min 0.836; 5 rows are closer to another head.

| class | rows | purity of the chunks (p) | identity to own head | reading |
|---|---|---|---|---|
| B2_Mm1a | 4 | 0.78–1.00 | 0.990–1.000 | clean |
| B2_Mm2 | 14 | 0.44–1.00 (10 of 14 ≥ 0.76) | 0.867–1.000; 1 of 14 closer to the Mm1t head | mostly clean |
| B2_Mm1t | 6 | 0.44–0.96 (5 of 6 ≤ 0.66) | 0.957–0.979; 3 of 6 closer to the Mm1a or Mm2 head | no chunk of its own: its copies sit between Mm1a and Mm2 |
| B3 | 27 | 0.56–1.00 | 0.914–0.990 | clean blocks, old (the copies are 23–27 % diverged) |
| B3A | 9 | 0.50–0.92 | 0.836–0.927; 1 of 9 closer to the B3 head | old and mixed with B3 |

The chunk order puts B3/B3A first (chunks 1–21), then B2_Mm2 → Mm1t → Mm1a (22–45), then B3 again
(46–60): the k-mer tree separates the B2 and B3 lineages and, inside B2, runs from the older Mm2 to
the younger Mm1a, with Mm1t in between and never pure.

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
B2a 22 rows (purity 0.72–1.00, 20 of 22 ≥ 0.92), B2b 2 rows (0.48, 0.60), B2c 2 rows (0.92, 0.48),
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

## Baseline: identity-threshold clustering (`vsearch_sweep.py`, output `vsearch_sweep.txt`)

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
