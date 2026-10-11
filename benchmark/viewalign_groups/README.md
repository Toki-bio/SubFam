# ViewAlign's group finding on SubFam's chunk alignments (2026-10-11)

Question (author): ViewAlign has clustering and a diagnostic-mutation search; can they define the groups better than the k-mer order and a fixed chunk size?
Answer by measurement: they work as a second stage on the **chunk-consensus alignment**, which is the documented ViewAlign workflow, and they give a hierarchy of groups
with named diagnostic columns. They need a well-aligned, shared region: on the ragged whole-copy L1 alignment they find coverage patterns, not subfamilies.

Code used: `cluster.js` of [Toki-bio/MSA-viewer](https://github.com/Toki-bio/MSA-viewer) at commit `5d637fb7c2d3b5576599225c8215562293ea371c`, class `SINEClusterer`
(position-pattern clustering = "Cluster Now"; `characterizeGroups` = exclusive / enriched columns of a given partition), loaded in node the way that repository's
`tests/clustering/*.test.js` load it (`cluster_rows.js`). UI defaults of those tests: minSize 3, minPerfect 5, quality 80/70/60 %, 10 iterations, `clusterChunked`.
Additions of ours: `--recurse 2` (re-run inside each cluster, same columns), `--min-occupancy 0.5` (mask columns that are a gap in more than half of the rows, see below),
`--groups` (characterise a given partition). Scoring: `score_clusters.py` (labels are counted in COPIES, a row = the copies of its chunk). Everything: `run_viewalign_groups.sh`;
results: `results/viewalign_groups.txt` and the JSON files. The viewer was not opened; its UI "Cluster Now" runs this same class on the same alignments:
`examples/b2_mm39_chr19_n50_chunks.aln.fasta` and `examples/l1pa_chr4_last1000_n20_chunks.aln.fasta` are the two input alignments (SubFam's own alignment of the chunk consensuses).

## B2 (60 chunk consensuses, 223 columns)
Top level, identical at every setting tried (minPerfect 5, 3, 2; with or without trimming and column masking): **two clusters**.
- 24 rows (chunks b2_22–b2_45, one run in the SubFam order; 11 exclusive columns): all B2_Mm2, B2_Mm1t and B2_Mm1a chunks; by copies B2a 1,080 of 1,200 (90 %), B2b 101.
- 36 rows (16 exclusive columns): B3, B3A, B2d, B2c, B2e chunks (B2d 1,167, B2e 409, B2c 195 copies).
That is the split of the 2021 scheme into the B2a/B2b lineage and the B2c/d/e lineage, and it puts Mm1a, Mm1t and Mm2 in one group as the 2021 paper does, with no label used.
Second level (minPerfect 2): the 24-row group splits into a Mm1a-like block (6 rows b2_34–b2_39, 2 exclusive columns; all 300 copies B2a), an Mm2-like block of 12 rows
(5 exclusive columns), inside which 3 rows (b2_43–b2_45, 2 exclusive columns) are the B2b-rich chunks (B2b 79, B2a 53 copies; purity only 0.53). The 36-row group **never splits**: 36 rows stay unassigned
at every setting, no exclusive pattern separates B3, B2d, B2e, B2c chunks. So an independent method with its own criterion (columns that belong to one group only) reaches the conclusion of the
chunk-purity analysis: the young indel-defined subfamilies have diagnostic columns, the old B2d/B2e do not. B2c (2 chunks) is below the 3-row minimum and is not found.
Given partitions (`characterizeGroups`, chunk-majority label): 2021 labels B2a 9 exclusive columns, B2d 21 (B2b, B2c: below 3 rows); RepeatMasker labels B2_Mm1a 0 exclusive columns, B2_Mm1t 0,
B2_Mm2 1, B3 10, B3A 4. Mm1a and Mm1t have enriched but no exclusive columns, in line with 2 nt between their consensuses.

## L1PA2–8
- Whole-copy chunk alignment (111 rows, 6,726 columns, `-c`): no useful result. The exclusive columns of small clusters are in the 5′ part that only a few consensuses cover (clusters of 3 to 5 rows with 248 or 283 "exclusive" columns); purity of the top-level clusters 0.27–0.66. Same cause as the length effect (`benchmark/line_real/README.md` §4–5): the alignment mostly records which part each consensus covers.
- Chunk alignment of the **last 1,000 bp** (111 rows; the alignment is gappy, 3,684 columns of which 160 are ≥ 90 % occupied): without a mask the clusters of 3 rows with 100–290 exclusive columns again come from sparse columns. With `--min-occupancy 0.5`
  there are 4 clusters and 17 unassigned rows: PA6 (6 rows, 0.92 of copies), PA7 + PA8 (24 rows, split into PA7 19 rows 0.90 and PA8 5 rows 0.72), PA3 (12 rows, 0.92) and a PA2/PA3 mix (15 rows, 0.50), PA4 (17 rows, 0.91) and PA5/PA6 (17 rows, 0.77);
  3 to 9 exclusive columns per group. This agrees with the 3′-end control (chunk purity 0.79) and gives named groups at PA3, PA4, PA6 and PA7 level; PA2 against PA3 and PA5 against PA6 do not separate by exclusive columns.

## Reading and cautions
- Use: SubFam (order, chunks, consensuses, alignment) then ViewAlign's clustering on the chunk alignment gives a short hierarchy of groups, each with its exclusive columns, instead of a fixed N; a group is "supported" when it has exclusive columns, and the 36-row B3/B2d group is the example of one that has none.
- It is not independent of SubFam's alignment: with a gappy alignment (L1), exclusive columns are artefacts of sparse columns; the occupancy mask is our workaround and is not a ViewAlign option (checked in `cluster.js` only; the UI may offer trimming).
- Labels (RepeatMasker, 2021 best hit) are used only to score, not to cluster. Rows are chunks (20–50 copies), so a mixed chunk puts all its copies into one group.
- The ViewAlign clustering is the author's own method (SINEClusterer); how it is named and cited in the paper is the author's decision (its features inventory calls it "the only comparable tool is command-line SubFam"). It also appears to share its logic with the author's unpublished peel / SINEderella line of work: the earlier decision not to expose SINEderella needs to be revisited before the paper uses this section.
