# Notes: Kryatova 2017, Konkel 2015, Gardner 2017 (MELT)

All three read in full from PDFs supplied by the owner (13, 15, 15 pp). Supplements (tables, alignments, supplementary methods) NOT read.

## Kryatova et al. 2017, Mobile DNA 8:6 - polymorphic AluS
- 52 reported polymorphic AluS elements examined; 48 confirmed AluS "using high stringency criteria".
- 49 in hg19 compared across RepeatMasker, the UCSC track and the 1000G (MELT/AluScan) call: 22 complete agreement, 22 minor disagreement (among AluS subfamilies), 5 substantial (1 AluSz vs AluJb, which became AluJ and was excluded; 4 AluS vs AluY).
- Resolution: five diagnostic nucleotides separate the six AluS subfamilies (Sc, Sg, Sp, Sq, Sx, Sz) from six AluY subfamilies (positions 95/93, 100/98, 197/196, 200/199, 219/218 in AluS/AluY numbering). Three elements were too 5'-truncated to be confirmed (two lacked the diagnostic positions; one carried a single diagnostic position matching AluY). One full-length element (11q14.1) carried AluS bases at six of twelve positions and AluY at three, others uninformative; classified AluS "due to predominating AluS features" (not a clean call).
- 37/48 (77 %) look like deletion polymorphisms, 7/48 (15 %) classical TPRT insertions (3' intact, poly-A, TSDs, L1 endonuclease site), 1 internal priming, 2 non-classical. These seven belong to four AluS subfamilies (Sx x3, Sg x2, Sp, Sz) and are significantly less diverged from their consensuses than all AluS in hg19 (mean 8.9 % vs ~11 %; permutation p = 0.0038).
- Relevance: subfamily labels of old (AluS) copies disagree between tools in more than half of the polymorphic cases tested; the resolution needed manual diagnostic-position reading, and truncation blocks classification.

## Konkel et al. 2015, GBE 7:2608 - 1000G pilot Alu insertions
- 343 (322 unique) Sanger-sequenced polymorphic Alu insertions, all AluY lineage, no AluJ/AluS. RepeatMasker: Ya5 48 %, Yb8 24 %, AluY 14 %, rest small. 96 % of copies <=2 % diverged from their consensus.
- High-resolution alignment against consensuses gave substructure ("bush-like", at least 42 active subfamilies; at least 28 potential source elements): 42 full-length AluY copies spread over about a dozen named lineages; Ya5 copies in 15 groups (Table 4: Ya5 103 of 149 full-length; Ya5a1 3, Ya5a2 6, Ya5b1 2, Ya5c1 4, Ya3.x, Ya4.x ...).
- Three newly named subfamilies, each from a single diagnostic change relative to a known one: Yb7a3 (3 sequenced copies; 8 exact hg19 matches), Yb8b1 (10 sequenced; 25 exact hg19 matches), Ya4a1 (11 sequenced; 13 exact hg19 matches; 24 independent insertions in total, which the authors say justify naming despite the variant being a CpG change).
- Their own words on nomenclature: more than one consensus has been reported under one subfamily name (Ya1, Yb10) or one consensus under different names (Yc, Yd).
- Found Yb8a1, Yb10, Yb11 (Ahmed 2013) once each in their data (locus 58, 325, 613).
- RepeatMasker label of locus 143 (all five Ya5 diagnostic changes) was not confirmed: "unable to identify a known consensus sequence available for this subfamily".
- 14 % (45/322) of copies 5'-truncated; two Yb8 copies too truncated to be authenticated.
- Supplementary files S1-S8 (alignments, GenBank KT305395-KT305737) not read. GenBank accession range is printed on p. 1.

## Gardner et al. 2017, Genome Res 27:1916 - MELT
- A mobile-element-insertion caller for short-read WGS, not a subfamily method; subfamily work is a side module.
- Subfamily module: CAlu finds interior mutations in Alu MEIs and assigns them to known subfamilies; LINEu does the same for L1. 79 "novel" Alu subfamilies were defined by shared interior mutations, requiring "at least five independent copies with a novel set of interior mutations (excluding CpG sites)" and a conservative threshold to remove Illumina errors; placed on a tree of known AluY subfamilies (Fig. 2). Details (supplementary methods, Table S8, Fig. S9-S10) not read.
- Interior mutations come only from reads supporting each insertion, so copies are partial and the definition is a count rule, not a statistical test.
- Population stratification: some subfamilies shared across all continental groups, others restricted (e.g. family F mostly African).

## What this means for us
- Orthology (human vs chimp flanks, polymorphism) is a separate dimension from the copy sequences SINEderella works with; the SINE-discriminator repo lists SINE_orth_loc / ComPair as the owner's own tools for it.
- Three independent groups define young subfamilies by one or two diagnostic changes plus a minimum number of independent copies (5 in MELT, 11-24 in Konkel for Ya4a1). A subfamily with that few copies is below SubFam's default chunk (50).
- Tool disagreement on old Alu labels is documented here (Kryatova) and in Carey et al.
