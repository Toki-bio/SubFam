# SubFam paper: draft abstract and claims (2026-10-08)

Status key: VERIFIED = measured and checked in this repository, PENDING = needs a run, NOT RE-RUN = reported in the repository, not repeated when this draft was written.
Numbers in [brackets] are placeholders for results that do not exist yet.

## Title (options)
- SubFam: bringing thousands of repeat copies to order and to the eye
- SubFam: threshold-free reduction of repeat families to an ordered alignment of consensuses

## Abstract (draft, ~190 words)
Repeat families, satellite arrays and multigene families can contain tens of thousands to millions of related copies, too many to inspect and too diverged for any single identity threshold. SubFam orders the copies along a k-mer guide tree, cuts the order into chunks of N neighbours, calls a plurality consensus for each chunk, and aligns the consensuses, producing a short alignment that can be read by eye while a table records which copies built each row. It needs no reference and no identity threshold; the only parameter, N, sets the resolution in copies. In simple simulations of eight subfamilies at 3, 8 and 15 % divergence (eight replicates each), SubFam recovered all eight source consensuses in every run with 100 rows for 2,000 copies; in a harder simulation with CpG decay, source elements, mixed ages and truncated copies it recovered 6.3-7.8 of 8 with 100 rows, at least as many as the best fixed VSEARCH identity in every scenario with 3 to 17 times fewer rows. The best VSEARCH identity moved with family age and the best identity threshold for VSEARCH moved with family age (0.95, 0.85, 0.70) and, in the oldest family, recovered 5.25 of 8 sources on average even when chosen per replicate using the truth. Rows are faithful to their copies in young and middle families and less so in old ones, where chunks mix sister subfamilies. A run on 100,000 copies took 30 minutes on four cores. On [N] hg38 Alu copies, the rows contain [x of y] Dfam and [x of y] published Alu consensuses at >= 98 % identity. SubFam does not call subfamilies or place individual orphan copies; it makes large repeat collections manageable.

## Claims (each with its evidence and what a reviewer will ask)

**C1. No identity threshold is needed to reduce a family of any age to rows that contain its source consensuses.** VERIFIED (simulation).
Evidence: benchmark/multiseed (seeds 1-8, 3 scenarios): `-n 20` recovered 8/8 in 24/24 runs with 100 rows; `-n 50` averaged 8.0 / 7.88 / 7.38 (young / middle / old) with 40 rows. VSEARCH --consout: the best single identity over the replicates was 0.95 (young, 8.0 recovered), 0.85 (middle, 8.0, 1,119 rows) and 0.70 (old, 5.25, 538 rows); with the identity chosen per replicate using the truth, old still averaged 5.25 of 8 (median 529 rows).
Reviewer: "VSEARCH recovers them too if you choose the identity." Answer: only inside an age-specific window (young 0.90-0.95, middle 0.80-0.85, old 0.70), with 92 to 530 rows; a mixed-age dataset has no single value. Do not say "parameter-free": N must be chosen.

**C1b. The advantage survives a harder simulation, at a lower level.** VERIFIED (4 seeds, benchmark/hard/README.md).
SubFam -n 20, CpG-aware identity: young 7.75, middle 7.5, old 6.25, mixed ages 7.0, 80 % truncated 6.5 (7.25 with -c), 100 rows. VSEARCH at its best fixed identity: 7.75 / 7.5 / 5.5 / 6.0 / 4.75 with 50 to 1,000 rows. Identity to the master over all columns is below 99 % for nearly every consensus in every method because CpG decay moves the consensus (state both). -c matters only when truncation is heavy; MAFFT guide-tree ordering (-m) fails under truncation (0.75 of 8), the k-mer tree does not.

**C2. N is an interpretable resolution parameter.** VERIFIED. A subfamily needs roughly N copies in the input to get a row of its own: `-n 50` loses the 40- and 60-copy subfamilies in some middle and old runs; `-n 20` does not. Cost: more rows.

**C3. The rows stand for their copies, and the paper says how well.** VERIFIED (simulation, seeds 1-4).
Evidence: benchmark/multiseed/fidelity.py: a copy's own row is its best-matching row for 98.1 / 87.3 / 57.5 % of copies (young / middle / old, `-n 50`; ceiling ~100 %), median identity to the own row 96.7 / 91.0 / 82.6 %.
Limit to state: fidelity falls with age because chunks mix sister subfamilies (chunk purity 0.99 / 0.91 / 0.71). Rows still contain their sources; the claim is organisation, not per-copy attribution.

**C4. The method scales to 10^5 copies on a workstation.** VERIFIED earlier (docs/SCALING.md), NOT RE-RUN here.
Evidence: 100,000 copies, `-P -n 50`, four cores: ordering 6 min 14 s, chunk alignments 6 min 12 s, final L-INS-i of 2,000 consensuses 17 min 35 s, total 30 min; 8/8 sources, purity 0.944. State the walls: the k-mer tree is n^2 (use it below about 20,000 copies, `-P` above); the final alignment grows as m^2 L^2 in the number of rows (hierarchical SubFam, i.e. SubFam on the consensuses, is the proposed route beyond that and is not yet measured).

**C5. The tool is a faithful, portable reimplementation.** NOT RE-RUN here (tests/check_equivalence.sh needs EMBOSS and node).
Evidence reported in the repository: k-mer order identical to ViewAlign's kmer-tree.js; the awk consensus byte-identical to EMBOSS `cons` over 3,600 runs; the MSF writer byte-identical to `seqret`. EMBOSS is not required.

**C6. On real data the rows contain the known consensuses.** PENDING.
Plan: benchmark/alu_hg38/ (protocol fixed in advance): recovery of Dfam Alu families and Price et al.'s 213 consensuses from three 30,000-copy samples, stratified by expected copies in the sample, against VSEARCH and `-m` baselines.
Illustration available now (not a validation): 316 Konkel Sanger Alu loci, 15 rows contain AluY / AluYa5 / AluYb8 at 1.000, 1.000, 1.000 identity, where VSEARCH needs about 140 rows and CD-HIT-EST about 125 for 0.99; the references are the consensuses used to label these copies (benchmark/alu_konkel/README.md).
Optional second real set: Timema SINEs with curated subfamilies (chunk purity 0.953 for 8 groups in the SINE-discriminator notes; labels partly derive from SINEderella's own assignment).

## Limitations to state in the paper
- Fixed-size chunks do not follow subfamily boundaries: redundant rows (an abundant subfamily fills several chunks) and mixed rows.
- Fidelity falls with family age (C3). Subfamilies smaller than N may be absorbed by a neighbour.
- The simple simulation is easy (independent decay from a master, no CpG, no truncation). The harder one (benchmark/hard) adds CpG decay, source elements, age mixture and truncation; there SubFam recovers 6.3-7.8 of 8, not 8/8. Both use parameters chosen by us, not calibrated to a genome, one family length (300 bp), no recombination.
- Truncated copies: consensus calling with the default plurality cuts a chunk's consensus to the region most copies cover; `-c` (coverage-relative plurality) exists but has only a small simulated check so far (a proper benchmark, e.g. the LINE simulation of docs/SCALING.md section 4, or leave it out).
- Time: 7-17 s for 2,000 x 300 bp on four cores, against about 1 s for VSEARCH.

## Not claimed
Subfamily identification or definition; placement or resolution of individual orphan copies; superiority to COSEG or to annotation pipelines (different tasks; see the SINEderella paper); that any reference annotation of Alu is right; uniqueness or priority.

## Related work to cite (read in this project)
Price et al. 2004 and COSEG (subfamily detection from alignments to one reference); Storer et al. 2021 (Dfam curation guidelines: COSEG, cd-hit); VSEARCH, CD-HIT, MMseqs2 (compared); ViewAlign kmer-tree.js (the ordering); MAFFT; EMBOSS cons.
