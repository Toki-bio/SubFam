# Mouse B2 SINEs (mm39 chr19): SubFam against two published classifications

Run 2026-10-10 in the cloud session (8 threads, MAFFT 7.505, SubFam 1.3.0). Everything is made by
`make_b2.sh WORKDIR` (downloads, extraction, sampling with seed 1, SubFam, grouped alignment).
Result file: `examples/b2_mm39_chr19_n50.aln.fasta`.

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
  continuum", and B3 as B2d. The B2a–e consensuses are not deposited; they are not in this file yet.

## Run
`SubFam.sh -n 50 -t 8 -x b2 b2_chr19_3000.fa`: 60 chunk consensuses, 19 s.
`build_examples.py ... --refs dfam_B2.fa --both`: for each rmsk class the consensus of all its
member copies (`TRUE_<class>_consensus_of_<class>_copies`, coverage-relative plurality 0.36), then
the Dfam consensus, then the SubFam rows whose majority is that class, in chunk order.

## What it shows (`rows.tsv`, from `examples/check_rows.py`)
Rows per class (majority): B2_Mm1a 4, B2_Mm1t 6, B2_Mm2 14, B3 27, B3A 9. Identity of the 60 rows
to the head of their own group: median 0.956, min 0.842; 6 rows are closer to another head.

| class | rows | purity of the chunks (p) | identity to own head | reading |
|---|---|---|---|---|
| B2_Mm1a | 4 | 0.78–1.00 | 0.990–1.000 | clean |
| B2_Mm2 | 14 | 0.44–1.00 (10 of 14 ≥ 0.76) | 0.867–1.000 | mostly clean |
| B2_Mm1t | 6 | 0.44–0.96 (5 of 6 ≤ 0.66) | 0.957–0.979; 3 of 6 closer to the Mm1a or Mm2 head | no chunk of its own: its copies sit between Mm1a and Mm2 |
| B3 | 27 | 0.56–1.00 | 0.914–0.990 | clean blocks, old (the copies are 23–27 % diverged) |
| B3A | 9 | 0.50–0.92 | 0.842–0.927; 2 of 9 closer to the B3 head | old and mixed with B3 |

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
  RepeatMasker; without the 2021 consensuses they cannot be marked in this file.
- Purity is against rmsk labels, which are themselves uncertain between Mm1a/Mm1t/Mm2 (2 and 16
  nt apart over 193 bp).
