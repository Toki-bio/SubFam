# Accuracy audit of MANUSCRIPT_BMC_v8.md, 2026-10-31

Four read-only audits (Abstract + Implementation; Examples 1-4 + legends; Tables 1-2 + simulation; Discussion + references),
plus DeepSeek's review. Every A-item was re-checked against the source file named. Line numbers refer to v8 at commit 5800f08.
Nothing below has been applied to the manuscript yet.

## A. Contradicted by the repo's own data (must fix)

| # | Line | Text now | What the data say | Proposed |
|---|---|---|---|---|
| A1 | 102 | "removing any step of the method loses the result" | ablation.tsv: no flip 0.771, MAFFT order 0.734 vs 0.793; only no-ordering (0.43) and no-chunks (0.42) lose it | "removing the ordering or the chunking brings chunk purity down to that of a single family consensus (0.42-0.43 against 0.79)" |
| A2 | 102 | "either fragments old families or merges young groups" | L85/Table 2: on the young Alu loci VSEARCH 0.95 matches SubFam (0.94 vs 0.92) | "either fragments the old families into singletons or merges across their groups; on a young family a working cutoff exists" |
| A3 | 45 | "label-free chunk statistics rank the variants in the same order as purity (Table 1)" | ablation.tsv: one consensus has mean_agree 0.796 > shuffled 0.757; Table 1 shows no such statistics | drop, or: "the lowest member identity (label-free) separates the ordered variants (0.71-0.72) from the unordered (0.59-0.60) and one consensus (0.50)" |
| A4 | 88 | "CpG decay keeps nearly every consensus of every method below 99 % identity" | hard/results.tsv: SubFam -n 20 recovers 7.75 (young), 4.75 (middle), 4.0 (mixed) of 8 over all columns; only the old family gives 0 | "CpG decay pulls the consensuses away from the master as divergence grows (old family: no source recovered by any method)" |
| A5 | 24, 99 | "4·n² bytes (0.8 GB at 10,000, 2.2 GB at 20,000, 3.6 GB at 30,000)" | 4·n² = 0.4 / 1.6 / 3.6 GB; 0.8 and 2.2 are measured peaks (scale_10k_20k.tsv 839, 2216 MB); 30k peaks 4.1-4.3 GB | "4·n² bytes (0.4, 1.6, 3.6 GB at 10,000, 20,000, 30,000 copies; measured peak of the whole run 0.8, 2.2, 4.1-4.3 GB)" |
| A6 | 31 | "B2_Mm1t ... has no rows of its own" | rows.tsv: Mm1t is the majority of 6 rows (b2_30, 33, 34, 39, 40 p0.96, 41); 3 of them closer to Mm1a/Mm2 heads | "B2_Mm1t is the majority of 6 rows, 5 of them at purity ≤ 0.66, interleaved with Mm1a and Mm2 rows" |
| A7 | 31 | "differs from B2_Mm1a at 2 of 193 positions" | Dfam rows in the chr19 alignment differ at 4 columns (29, 135, 147 substitutions; 64 N in Mm1t) | "at 3 of 193 positions (and one N)"; re-check against Dfam |
| A8 | 39 | "A genome-wide Alu run, on 30,000 copies of hg38 chr1" | alu_hg38/README: one chromosome, one sample; genome-wide run was left for the servers | "A run on 30,000 Alu copies of hg38 chr1" |
| A9 | 42 | "in 6.3 min at 4.1 GB ... (13 min with the processor shared)" for seed 1 | runs.tsv: seed 1 timed once, 791 s shared CPU; 6.3 min is seeds 2-3. -P (439 s) and Alu times also shared CPU, so the text makes -P look slower | "13 min at 4.1 GB on a shared processor (6.3 min on an idle one for two further samples) and 7 min at 0.7 GB with -P (shared)", or re-time on an idle machine |
| A10 | 22 | "The default of 50 is the smallest chunk that still keeps the properties of a sequence group" | sensitivity.tsv: N 20/30/50/100 purity 0.812/0.803/0.793/0.757; no test finds a minimum | "N = 50 by default: at p = 0.36 a base then needs 18 agreeing copies; smaller N gives purer but more numerous and noisier rows" |
| A11 | 22 | "rows of N = 50 represent their copies more closely than those of N = 20" | fidelity.txt: middle family N = 20 is better (88.0 vs 87.3) | add "equally in the middle one (87.3 against 88.0 %)" |
| A12 | 99 | "each real example is one chromosome of one genome" (also L85 "each family is one chromosome") | Example 4 is chr1-3; Example 3 is genome-wide loci | "each real example comes from one genome" |
| A13 | 95 | "a different k for old ones" | SubFam.sh:48 "best at all divergences tested"; README k = 5-6 best everywhere | delete |
| A14 | 26, 161 | "Every example alignment is one MAFFT alignment of three kinds of row" | 30k B2 alignment has no heads and no rows.tsv; Konkel has no CLASSCONS; L1 CLASSCONS from ≤ 60 random members | "Most example alignments ..." + state the exceptions |
| A15 | 110 | "MAFFT ≥ 7.5" | SubFam.sh:98, README: mafft ≥ 7; 7.505 is the tested version | "MAFFT ≥ 7 (7.505 tested)" |
| A16 | 36 | "17 of 111 rows lie closer to a neighbouring class's head" | rows.tsv: l1_110 PA8 is closer to PA6, two classes away | "closer to another class's head, 16 of them a neighbouring one" |
| A17 | 80 | Table 2, Alu vsearch 0.99: 35 % | 109/316 = 34.5 % → 34 %; recompute 96.5/28.5/15.5 % cells from integer counts | fix rounding |

## B. Overstated or not shown (should fix)

| # | Line | Text now | Problem | Proposed |
|---|---|---|---|---|
| B1 | 11, 85 | "fragments the old families into singletons or merges across their groups" | b2_sweep.txt: vsearch 0.60 is intermediate (63 % covered); **under the Vassetzky 2021 labels its groups are purer than SubFam's (0.849 vs 0.773)**, 0.70 gives 0.954 at 40 % | "no cutoff gives groups that are both pure and cover all copies"; report the 2021-label purity for B2 |
| B2 | 11 | "on the young family it reaches the same purity ... while leaving a tenth of the copies unplaced" | the 11 % are in groups of 2-9; 8 singletons (2.5 %); only VSEARCH matches, MeShClust2 does not | "on the young loci VSEARCH reaches SubFam's purity with a tenth of the copies in groups of fewer than ten" |
| B3 | 11 | "N ... is its one substantive parameter" | k moves purity by 0.04, N by 0.05; -c, -r matter for their inputs | "its main parameter is N" |
| B4 | 25 | "8,000 random alignments under two awk implementations ... chunk splitter ... 1 to 30 and N 2 to 10 (code review of 2026-10-11)" | not found in the repo (only review is CODE_REVIEW.md 2026-10-10, v1.2.0; the 8000 in tests/ is the k-mer case generator). What exists: check_equivalence.sh, 3,600 consensus runs, 80 MSF, 0 differences | add the tests to tests/, or cite what exists. **Author: if these results are elsewhere, say where** |
| B5 | 22 | "N = 100 dilutes the rare groups further" | rare groups measured only at N 20 and 50 | "N = 100 lowers chunk purity on B2 (0.757 against 0.793)" |
| B6 | 21 | "(most of the gain at a tenth of the time)" | only a code comment; no timing of the two-round mode | delete or give the 0.004 identity measurement |
| B7 | 23 | "600 rows is about the most ... can take in"; "close groups compress to the same few hundred consensuses" | no test; 588 of 600 rows distinct at 30k | give 600 as the author's practical limit; drop "compress" |
| B8 | 42 | "so the picture of Example 1 holds at ten times the size" | measured: equal purity and same block layout; only 11 of 36 chr19 B3/B3A rows recovered | "chunk purity is the same as at 3,000 copies (0.805 against 0.793)" |
| B9 | 42 | "35 of the 60 rows ... the rest being B3/B3A rows" (DeepSeek misread it as a contradiction; it is not) | recomputed: all 24 B2 rows + 11 B3 rows; no B3A row; other samples 37-40 | "all 24 B2 rows (chunks 22-45) and 11 of the 27 B3 rows ...; none of the 9 B3A rows (35 of 60; 37-40 in other samples)" |
| B10 | 42 | "the structure of the 60-row run at ten times the resolution" | same N, so same resolution | "in ten times as many rows" |
| B11 | 36 | "chunks follow them at the resolution of neighbouring classes (PA2 + PA3, PA4-PA6, PA7 + PA8)" | grouping measured on peel groups (SubFam 1.2.0 + peel), not on chunks | "the chunks mix neighbouring classes (92.6 % of copies within one class of their chunk's majority)" |
| B12 | 31 | "B3A ... shares its rows with B3" | B3A is the majority of 9 rows, purity 0.50-0.92 | "is the majority of 9 rows, mostly impure, interleaved with B3 rows" |
| B13 | 33, 42 | "about 70 % identical"/"25-30 % divergence between copies"/"copies are 25-30 % diverged" | measured: 30-33 % from consensus (Vassetzky labels); no pairwise value | "about 30 % diverged from their consensus" in both places |
| B14 | 33 | "this is the limit of every copy-level method" | only three methods tested | "neither SubFam nor identity clustering (Table 2) separates them" |
| B15 | 91 | "The ablation locates the contribution ... in the fixed chunk" | no variable-chunk variant in Table 1 | "in the combination of ordering and chunking" |
| B16 | 102 | "the rows carry the known consensuses of the family" | B3 0.943, B3A 0.889; Mm1t; L1 17/111 | "carry the consensuses of the young classes (0.98-1.00) and approach those of the old ones (0.89-0.94)" |
| B17 | 102 | "the input that curation ... need and cannot take from the raw copies" | not shown; "subfamily definition" twice | "a compact input for curation, subfamily definition and phylogenetic analysis, which are hard to apply to thousands of raw copies" |
| B18 | 97 | "SubFam serves as a step of input preparation for them" (incl. COSEG, Price) | COSEG/Price work on copies, not consensuses; not shown | "Each answers a different question; SubFam's consensuses are a compact input for downstream curation and alignment" |
| B19 | 95 | cascade: "-e gives the copies ...; the remainder is the next input ... every copy accounted for" | -e writes no remainder; the procedure was not benchmarked (deplete_loop.sh is a different, automatic one) | "(the remainder, built with any FASTA tool, ...)"; say the cascade is proposed, or cite SCALING.md §6 result |
| B20 | 24, 108 | "stops with a message when it does not fit"; "Linux, macOS" | check reads /proc/meminfo (Linux only); real peak 1.2-2x the matrix; no macOS run recorded | "On Linux SubFam compares the matrix size with the free memory ..."; "Linux (tested); other Unix with bash should work" |
| B21 | 88 | "10 of the 40 rows ... an abundant subfamily filling several chunks" | one run, middle family, seed 1; six subfamilies repeat; not in a benchmark file | "In one run (middle family, seed 1) ... each abundant subfamily fills several chunks" |
| B22 | 88 | "VSEARCH identity that recovers the sources ... 0.70 for the old one" | at 0.70 it recovers 5.25 of 8; middle omitted | "best VSEARCH identity moved with age (0.90-0.95, 0.85, 0.70, where it recovered 5.25 of 8)" |
| B23 | 120 | "every example and benchmark is reproduced by ..." | alu_hg38/run_chr1_30k.sh, multiseed, run_scale.sh, sensitivity.sh, seeds.sh, hard/, alu_konkel not listed | list them, or "each benchmark/ folder has its script" |
| B24 | 158, 165-171 | Figure 2 "220 of 1,105 columns"; Additional files 2, 4, 6 | column count only printed to stderr; the files 2/4/6 are not built yet | commit the figure log; build the files |

## C. Wording and clarity (incl. DeepSeek points checked against the code)

- L21 chunk count: code forces at least one chunk (SubFam.sh:565) → "round(n / N), at least one".
- L21 "closest ends": define: the four distances between the end copies of the two subtrees; the pair with the smallest distance meets (SubFam.sh:351-356).
- L21 threshold: "positive-scoring matches" = sequences whose base scores positively against the chosen base under EDNAFULL (ambiguity codes count, gaps never).
- L21 with -c, lower case still uses half of all chunk members (SubFam.sh:446, 480): say so.
- L21 final alignment rows are in MAFFT --reorder order; the consensus file and chunk table keep k-mer order (SubFam.sh:702). L42 and Additional file 7 say "SubFam's row order".
- L21 "MAFFT L-INS-i ... and in MAFFT's automatic mode beyond": --auto is not L-INS-i.
- L21 chunk statistics: agreement/identity are against each column's most frequent base, not the consensus; ragged = share of core columns with ≥ 20 % gaps.
- L25 the 26 robustness cases include malformed inputs; refusal tests are check_inputs.sh.
- L23 my sentence: "seed 1" → "(seed 1; seeds 2 and 3 for the repeat samples)"; "examples use 3,000 and 30,000" omits Examples 2-3.
- L42 "Full-length B2 copies" → "B2-family copies of at least 170 bp".
- L42 B3 0.943 is seed 1 only; seeds 2-3 give 0.954/0.952.
- L52 "same chunks" → "same chunk size".
- L64/83 Table 2 legend: name the labels (RepeatMasker for B2/L1, Price for Alu).
- L85 "COSEG, run once" → "run on 295 of the loci at three minimum group sizes".
- L88 "98.1, 87.3 and 57.5 %" add "(N = 50; seeds 1-4)".
- L91 "built without a cutoff" → "without an identity cutoff".
- L97 COSEG: "copies that span a common window of one reference alignment".
- L99 "(0.472 against 0.400)" → "(0.400 against 0.472)"; "lifts" → "needs far less memory (0.8 GB at 30,000)".
- L99 "3.6 GB" → matrix 3.6 GB, run peak 4.1-4.3 GB.
- 3000/2222 in tables vs 3,000/2,222 in text; CpG and UCSC missing from abbreviations.
- Table 2 bold: present in markdown and carried into the docx (checked).

## D. Repo, not manuscript

- SubFam.sh:566-567 comment says 0.399 vs 0.466, 3.5 GB; manuscript/README 0.400/0.472, 4.3 GB.
- sine_b2/README: "no chunk of its own" for Mm1t contradicts its own "6 rows" column; "2 and 16 nt apart".
- alu_konkel/README reports 15 chunks (SubFam 1.2.0); the manuscript's 16 is current.
- REFERENCES (d)2 (why 316 of 322 Konkel loci) and (d)7 (Dfam 4.0 citation) still open.

## E. Waiting on the author

ORCID, Sergei/Sergey, affiliation, e-mail, Zenodo DOI, funding, competing interests (co-author of Vassetzky 2021), acknowledgements, LLM statement, prior-use papers (renumbers MAGUS), the 8,000-alignment test (B4).

## Status 2026-10-31

A, B (incl. B1 with the Vassetzky-label purity reported in the Table 2 text and in Limitations), C and D applied. B4: the unverified 8,000-alignment and splitter sentence replaced by check_equivalence.sh (3,600 / 80 runs) and check_inputs.sh. Still open: Figure 2 column count log (B24), Additional files 2, 4, 6 not built, REFERENCES (d)2 and (d)7, section E.
