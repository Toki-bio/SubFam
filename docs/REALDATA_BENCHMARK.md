# Real-data benchmark plan (draft)

Goal: show on real genomes that SubFam's chunk consensuses recover biologically meaningful
subfamilies, without trusting any single annotation as copy-level truth.

## 1. What the literature says about Alu labels

Source: a literature pass done from search-result text only. Several publisher sites
(Genome Research, PMC, Dfam) could not be fetched, so items marked *check* need verification
against the full paper before they go into a manuscript. No paper was found that quantifies
the copy-level error rate of Alu subfamily labels genome-wide; we have to measure it.

| Level | Basis | Trust |
|---|---|---|
| AluJ / AluS / AluY | Jurka & Smith 1988; many shared diagnostic substitutions | high (stable since 1988; error rate not measured, *check*) |
| Young Y lineages (Ya5, Yb8, Yb9, Yb10, Yb11) | diagnostic mutations plus source-copy tracking (Ahmed et al. 2013, Mobile DNA 4:25: one Yb8 copy sourced 712 copies, 54 % of human-specific Yb8) | high, with caveats below |
| Y lineages with several concurrent sources (Yh, Yi) | Styles & Brookfield 2009, BMC Evol Biol 9:102 | medium |
| MELT "novel" subfamilies (79 reported) | >= 5 independent copies with a novel non-CpG interior mutation set (Gardner et al. 2017) | medium; data-driven, not curated |
| AluS and AluJ sub-subfamilies (Sx, Sx1, Sx3, Sg, Sp, Sq, Sc, Jo, Jb) | consensus similarity of old, CpG-decayed copies | low |
| Any Dfam/RepBase entry without independent support | library curation unknown | do not use as truth |

Documented failure modes:

- RepeatMasker, UCSC and 1000 Genomes labels agreed fully for only 22 of 49 AluS polymorphisms; disagreements were mostly S vs Y, and three 5′-truncated copies could not be classified (PMC5402677).
- Some Ya5-lineage insertions carry only 2 or 3 of the 5 Ya5 diagnostic changes; one with all five was labelled AluYk3 by RepeatMasker (Genome Biol Evol 2015, 7:2608).
- Gene conversion can replace an Alu body by a younger subfamily's sequence while the flanks keep their old phylogeny; body-only calls are then wrong. Reported to contribute 10-20 % of variation among recent Ya5-related elements (*check* which paper).
- Price, Eskin & Pevzner 2004 (Genome Res 14:2245) found 213 statistically separable subfamilies from the 31 in Repbase at the time, and their consensuses differ from Repbase mainly at CpG sites: sub-subfamily structure is closer to a continuum than to discrete classes.
- Labels depend on the library version, so every label must be stored with the exact Dfam / RepeatMasker version.
- Not covered by the pass: Konkel 2010, Salem 2003, Kojima, any 2020-2025 T2T or pangenome re-derivation of Alu subfamilies, AluScan-type classifiers.

## 2. Consequences for the design

1. No flat truth. Score at three levels and report each separately:
   (a) J/S/Y (robust), (b) Y lineages (Ya5, Yb8, ...) from diagnostic positions and source tracking, (c) AluS sub-labels as "soft", never as pass/fail.
2. Use two independent label sources per copy: library label (RepeatMasker with a stated Dfam version) and a diagnostic-position call computed by us from the alignment. Copies where they disagree go to an "ambiguous" class and are excluded from purity, but their fraction is reported.
3. Exclude or flag: 5′-truncated copies (below a fraction of the consensus length), partial diagnostic matches, and body/flank discordance where an orthologous-flank age is available.
4. Report per-level: purity and recovery of each reference subfamily (a reference subfamily is recovered when some chunk consensus has it as plurality of its copies and is closer to its consensus than to any other), adjusted Rand index of chunk-vs-label, and the number of redundant and mixed rows. Compare to VSEARCH (best per-family `--id`), CD-HIT, MMseqs2, same inputs.
5. Add a simulation layer for exact truth: multi-source tree, CpG-biased mutation, back-mutation, gene conversion between lineages. Real data cannot give exact truth, simulation can (extends benchmark/simulate.py).

## 3. Datasets, in order of effort

| Case | Data | Truth | Purpose |
|---|---|---|---|
| A. Young Alu (Y lineages) | hg38 / T2T-CHM13 AluY copies >= 250 bp, plus 1000G polymorphic insertions | diagnostic positions + Ahmed 2013 Yb clusters; MELT subfamily field (*check VCF INFO*) | main Alu result |
| B. Old Alu (J, S) | genome-wide AluJ/AluS copies | J/S/Y only; sub-labels soft | scale and old-family behaviour |
| C. Young L1 (L1HS, L1PA2-L1PA8) | 5′UTR/3′UTR diagnostic changes; use the common 3′ span | subfamily labels (*check*) | cleaner small-scale truth, truncation test for `-c` |
| D. Satellite | alpha satellite HORs in T2T-CHM13 (monomers, HOR variants) | HOR structure annotation | non-TE case |
| E. rDNA / multigene family | rDNA units or a large gene family | known variant classes | multigene case |
| F. Simulation | benchmark/simulate.py extended | exact | sensitivity to CpG, back-mutation, conversion |

## 4. Open decisions

- Which Dfam version(s) to freeze (record version with every label).
- Whether to include L1/SVA in the first paper or keep Alu + one satellite + simulation.
- Whether a Price-2004-style tree on the same copies is reproducible as a second independent labelling (tool availability unverified).
- Verify every *check* item above against the full texts before citing.
