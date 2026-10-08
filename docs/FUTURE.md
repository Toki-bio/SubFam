# Future steps by element class, and FaSort10 as a pre-sorter (2026-10-08)

Rule for everything below: a step is "done" only when the claim has a test with a score fixed in advance, first on simulation, then on real data with a truth that does not come from the method.
Status words: DONE (in this repo), KNOWN (measured elsewhere in the repo, see the file named), IDEA (not measured). Statements about tools or datasets I have not opened are marked UNVERIFIED.

## 0. What the paper's claim needs from each class
SubFam brings an unorganised multitude to order and to the eye. For each class the question is the same: does the short alignment still hold the known consensuses, and is the provenance (which copy built which row) still usable?
What changes between classes is what makes copies comparable: shared length (SINE), truncation (LINE), two similar ends (LTR), tandem order (satellite), module content (mosaic).

## 1. LINEs
Known (docs/SCALING.md section 4, simulated 6 kb families): full-length copies order at chunk purity 0.99; the same copies 5'-truncated order at 0.45; Jaccard on a common 2.1 kb span lifts it to 0.77; sister subfamilies differ at 1.3 % of sites, so fragments of 300 bp carry about 4 informative sites and cannot be placed by any method.
Known (benchmark/hard, 80 % truncated 300 bp copies): `-c` lifts recovery 6.5 to 7.25 of 8; the MAFFT guide-tree order (`-m`) collapses (0.75); the k-mer tree does not.
Next steps:
1. **A LINE simulation that looks like L1**: 6 kb master, 5'-truncation with an exponential length distribution, a fraction of 5'-inverted copies (the MELT paper reports 18.2 % of L1 insertions with a 5' inversion, and that inversion breakpoints are common, so strands are not uniform inside one copy), a CpG-rich 5' UTR, internal deletions, subfamilies defined in the 5' UTR and the 3' end. Extend simulate2.py; score with evaluate2.py plus a coverage-per-row profile.
2. **Discover from long copies, assign fragments** (SCALING section 5, "SubFam-XL"): implement as a documented recipe first (length quantile, trim to common span, SubFam on the long copies, assign fragments to the rows), then measure how many fragments land in the right row and how many rows change when fragments are added.
3. **Position-aware k-mer distance** on the shared interval instead of whole-sequence Jaccard (SCALING section 4, design consequence 2). Test on the simulation above against the current ordering.
4. **Real data**: hg38 L1 copies with RepeatMasker/Dfam labels (L1HS, L1PA subfamilies; UNVERIFIED which label set), restricted to the 3' 1-2 kb where most copies are complete; metric: rows containing the Dfam consensus of that region at a fixed identity, as in benchmark/alu_hg38. Truth caveat: labels come from alignment score (Carey et al. 2020).
5. Scale: full-length L1 is about 6 kb; MAFFT per chunk is O(L^2) (fine at 50 x 6 kb, seconds). The 5 kb limit of VSEARCH's SIMD aligner matters only for the VSEARCH baseline.

## 2. LTR retrotransposons (LTR-internal-LTR, solo LTRs)
Known (Storer et al. 2021, read): Dfam's curators use cd-hit plus a script for solo LTRs and internal deletions because COSEG-type alignment structure does not fit them; the cd-hit example reached 94 of 133 solo LTRs in clusters.
Next steps (all IDEA):
1. **Treat the two parts separately**: extract LTR and internal sequences as two inputs (annotation tells where they are), run SubFam on each, keep the locus coordinates so a copy's LTR row and internal row can be joined afterwards. SubFam never sees both in one sequence.
2. **Simulation**: LTR 400-800 bp, internal 5-8 kb, solo-LTR share 80-90 %, internal deletions, both LTRs of an insertion identical at birth (age = their divergence), subfamilies defined in the LTR (as is common for HERVs; UNVERIFIED for specific families). Score LTR-row recovery and the fraction of solo LTRs placed.
3. **Dating as an extra output**: for paired LTRs the divergence between the two LTRs of one copy dates the insertion; SubFam could report, per row, the distribution of that divergence (needs the pairs; nothing to build in SubFam itself, an annotation step).
4. **Real data**: HERV families in hg38 from Dfam (UNVERIFIED which and how many); metric as for Alu.

## 3. Satellite monomers and their stacks (arrays, higher-order repeats)
Known: the repository's own notes (SINEderella docs/SATELLITES.md, docs/ARRAYS.md; not read in this session) deal with arrays and SINE-derived satellites inside SINEderella.
The new problem is order along the array: neighbouring monomers in an array are far more similar to each other than random copies are (higher-order repeat copies can be nearly identical), so a chunk made of nearest neighbours may be one stretch of one array, not a variant class.
Next steps (IDEA):
1. **Simulation**: monomers of about 170 bp with a few monomer types; higher-order repeats built by stacking k monomer types in a fixed order; arrays made of HOR copies with 1-3 % divergence between copies of the same array and a larger divergence between arrays; monomer-level mutation. Truth: monomer type, HOR unit, array.
2. **Sampling that respects arrays**: take at most m monomers per array (or cluster adjacent monomers first), as SINEderella's array flag does for SINEs; compare SubFam rows with and without it.
3. **Two-level use**: SubFam on monomers gives monomer rows; re-encode each array as a string over the row alphabet (row index per monomer); periodic structure of that string is the HOR; run a second, different reducer on those strings (not SubFam as it is: the objects are strings of symbols). Measure whether the HOR unit is recovered from the string.
4. **Scale**: alpha satellite alone runs to hundreds of thousands of monomers in a human assembly (order of magnitude, UNVERIFIED); that needs the `-P` route or partition-then-SubFam (SCALING section 3) and assignment of all monomers to rows afterwards.
5. **Real data**: T2T-CHM13 alpha-satellite monomer and HOR annotation (UNVERIFIED: where the annotation is and in what format; to be located). Metric: monomer types contained in the rows.

## 4. Mosaic and composite elements
Known (SINEderella docs, read): composites (a unit plus a second unit, such as r5_r6_P26 in the rsi bank) account for 83 % of vote failures in the rsi run (docs/FAMILY_SUBFAMILY_ASSIGNMENT.md); RESEARCH_DIRECTIONS.md proposes a module library plus a presence/absence matrix ("panconsensus").
What goes wrong in SubFam (from reading SubFam.sh): copies with different module content are aligned by MAFFT in one chunk, the plurality consensus is called column by column, and a chunk that mixes arrangements gives a consensus that does not exist as a real copy.
Next steps (IDEA):
1. **Detect, do not hide**: report per row a coverage profile (copies covering each column, which `-c` already computes internally). A row whose coverage is bimodal (a segment present in some copies and absent in others) is flagged "mixed module content". Cheap, and it turns the failure into information.
2. **Simulation**: modules A, B, C of 100-150 bp, composites AB, BC, ABC and singles, module-level divergence, 10-30 % of copies truncated at a module boundary. Score: recovery of the module consensuses, fraction of rows that are chimeric, fraction of mixed rows that the coverage flag catches.
3. **Module-level mode**: split copies into windows with a segmenter (SINEderella's flankscan/segmentation or breakpoint windows), run SubFam on the windows pooled, return each copy as a path through module rows (a presence matrix). This is the panconsensus idea with SubFam as the module-library builder; keep it as a separate workflow, not an option inside SubFam.
4. **Real data**: the rsi composites and the eri case in SINEderella's notes (bimodal flank population) are the only cases with an owner-made truth.

## 5. FaSort10 as a smart pre-sorter
What FaSort10 is (repo Toki-bio/FaSort10, 54 lines read in full, plus a copy of asSINEment): given a set of group consensuses ($1) and a bank of sequences ($2), it splits the bank into 20,000-sequence parts (seqkit), runs `ssearch36 -g -3 -z 11 -E 2 -w 95 -W 70 -m 8` of the consensuses against the bank ten times, keeps for every sequence its 10 best-scoring hit lines over all cycles, tallies which group each belongs to (output `seqid g1_7_g2_3`: group and number of its lines among the top 10), writes per-group id lists and per-group FASTA banks (`.bnk`). It is the ancestor of asSINEment, without asSINEment's unanimity rule and 0.45 threshold.
Practical notes from the code: threads are fixed at 24; part files are found with a `*.part_*` glob and scores are appended to `score` in the working directory (rerunning in the same directory accumulates); the shebang is `#/bin/sh` (a missing `!`) while the body uses `[[` and gawk-only arrays of arrays.
Where it fits SubFam (my reading of "smart sorter": sort the multitude into coarse bins using given prototypes before or after SubFam; correct me if you meant something else):
1. **Partition first at scale.** SCALING section 3 partitions with a loose VSEARCH pass; FaSort10-style Smith-Waterman voting against a handful of coarse prototypes (for instance the consensuses of a first SubFam run on a 10 % sample) is slower but does not lose copies below about 75 % identity, where VSEARCH's k-mer candidate ranking fails (SCALING section 2). Test: old and mixed scenarios of benchmark/hard, plain SubFam against sample -> prototypes -> sort all copies -> SubFam per bin; score recovery, rows, time.
2. **Assign every copy to a row afterwards.** Chunk membership is what `chunks.tsv` records; sorting all copies against the final rows gives a chunk-independent membership, the vote counts show which copies are ambiguous, and orphans become a visible list (the "poor orphan copies" of the paper's narrow claim). Test: fidelity.py on benchmark/hard before and after.
3. **Two-level pipeline for 10^6 copies.** Level 1: SubFam on a sample; level 2: sort all copies into the level-1 groups; level 3: SubFam inside each group; pool the consensuses. This is SubFam-XL (SCALING section 5) with FaSort10 as step 5.
4. **Fragments (LINE).** Assign short fragments to the rows built from long copies, as the LINE recipe in section 1 needs.
Implementation sketch: an optional `--presort GROUPS.fa` in SubFam.sh that calls a vote script (ssearch36 if installed, `vsearch --usearch_global` otherwise) and runs the main pipeline per bin; not before the experiment in item 1 shows a gain.

## 6. Order of work
1. Mosaic coverage flag (small, turns a failure into output) and the LINE simulation (section 1, step 1), because the paper's limitations section names truncation.
2. hg38 Alu consensus recovery (benchmark/alu_hg38), the paper's real-data claim.
3. FaSort10 experiments (section 5, items 1 and 2), then decide whether the pre-sort belongs in the tool or in SINEderella.
4. LTR and satellite simulations, then real data; mosaic module mode last.
