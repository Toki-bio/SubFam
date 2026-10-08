# Konkel et al. 2015 Alu loci: a first real-data check (coarse)

Data: the 343 GenBank records KT305395-KT305737 (Sanger-sequenced polymorphic Alu insertions, 1000 Genomes pilot; Konkel et al. 2015, GBE 7:2608). Download them yourself with efetch (see docs/HANDOFF_vscode.md); the records are not stored here.

`prep.py GB_DIR PRICE_SUPPLEMENT OUT_DIR`
1. local-aligns every locus to the Price et al. 2004 AluY consensus (supplement of that paper, not stored here), cuts the Alu body out;
2. labels each body with the best-scoring of four Price consensuses (AluY, AluYa5, AluYb8, AluYb9); a tie is "ambig";
3. writes `alu_bodies.fa` (316 full-length bodies, >=90 % of the AluY consensus, 250-330 bp) and `loci_labels.tsv` (all 343 loci, committed).

Result of `SubFam.sh -n N alu_bodies.fa` (k-mer ordering, defaults otherwise), chunk majority purity against those labels:

| -n | chunks | purity | purity of random chunks of the same sizes |
|---|---|---|---|
| 20 | 15 | 0.918 | 0.53 |
| 10 | 31 | 0.937 | 0.55 |

Labels: 174 AluYa5, 81 AluYb8, 73 AluY, 15 ambiguous among the 343 loci (Konkel et al. report Ya5 48 %, Yb8 24 %, AluY 14 % with RepeatMasker on 322 unique loci). Pure chunks cover almost all Ya5 and AluY copies; the mixed chunks sit at the AluY / Yb8 border (a chunk of 8 Yb8 + 7 AluY + 4 Ya5 at n=20).

## What this does and does not show
- It shows that k-mer ordered chunks recover the three coarse Y lineages on real Sanger loci, far above random chunking.
- It does not show subfamily recovery at Konkel's fine level. The labels come from alignment score to consensuses, the same kind of similarity the ordering uses, so the two are not independent. Truth at the level of Yb7a3 (3 copies), Yb8b1 (10), Ya4a1 (11), Yb10 (8) and Yb11 (16) needs the authors' per-locus calls (their supplementary file S3, not available to us) or diagnostic-position calls made by us.
- Each locus is a Sanger consensus of one individual's allele; 343 records contain 322 unique loci (duplicates between pilot sets).
- No separation claim for AluS or AluJ: these data contain none.

## COSEG on the same loci (single run, 2026-10-08)
COSEG (github.com/rmhubley/coseg, Price's method, `runcoseg.pl -k -d -m M`) needs copies aligned to one reference. `to_coseg.py` aligned the 316 bodies to the AluSx consensus shipped with COSEG and dropped copies missing >5 bases at an end (Price's rule): 295 copies remain (our conversion, not COSEG's own preprocessAlignments.pl). Same 295 copies for SubFam. `compare.py` scores both against the lineage labels (11 ambiguous copies removed):

| method | groups | purity | homogeneity | completeness | V |
|---|---|---|---|---|---|
| COSEG -m 50 | 2 | 0.778 | 0.522 | 0.966 | 0.678 |
| COSEG -m 10 | 3 | 0.778 | 0.522 | 0.776 | 0.624 |
| COSEG -m 5 | 4 | 0.806 | 0.559 | 0.705 | 0.623 |
| SubFam -n 100 | 2 | 0.729 | 0.366 | 0.588 | 0.451 |
| SubFam -n 75 | 3 | 0.768 | 0.507 | 0.493 | 0.500 |
| SubFam -n 50 | 5 | 0.778 | 0.585 | 0.377 | 0.458 |
| SubFam -n 20 | 14 | 0.972 | 0.928 | 0.356 | 0.515 |

Reading (what the numbers support and what they do not):
- COSEG found a Yb8 group against the rest (66 Yb8 together; 10 ambiguous copies went with it) but never separated Ya5 from AluY (217 copies together) even at -m 5; the Ya5 pattern (five diagnostic changes) was not split off here. Reason not investigated.
- SubFam chunks of fixed size give homogeneous but numerous rows (homogeneity 0.93 at -n 20, completeness 0.36): by design one lineage fills several chunks. At the same granularity (2-5 groups) COSEG scores higher than SubFam chunks; at -n 20 SubFam is far more homogeneous but over-split.
- SubFam chunks are not the end product of the SINEderella workflow (a person or the peel merges them), so this compares COSEG with the chunk step only, not with the peel.
- Limits: 295 copies, three coarse labels derived by similarity (not independent truth), one COSEG run per setting, minimum subfamily size below COSEG's design point (-m 50 was meant for ~480,000 copies), our own alignment conversion.

## Worked example: rows needed to contain a known consensus (2026-10-08)
Same 316 Konkel Alu bodies. For each method, the identity of the best output sequence to four Price et al. consensuses
(`recover.py`, coverage >= 0.90). SubFam rows are consensuses; VSEARCH and CD-HIT-EST rows are single copies (centroids / representatives).

| method | output rows | AluY | AluYa5 | AluYb8 | AluYb9 |
|---|---|---|---|---|---|
| SubFam -n 20 | 15 | 1.000 | 1.000 | 1.000 | 0.993 |
| SubFam -n 10 | 31 | 0.997 | 1.000 | 1.000 | 0.993 |
| VSEARCH --id 0.90 | 4 | 0.965 | 0.982 | 0.976 | 0.969 |
| VSEARCH --id 0.95 | 18 | 0.975 | 0.982 | 0.976 | 0.969 |
| VSEARCH --id 0.98 | 70 | 0.993 | 0.989 | 0.990 | 0.990 |
| VSEARCH --id 0.99 | 140 | 0.996 | 0.996 | 0.997 | 0.993 |
| CD-HIT-EST 0.90 / 0.95 / 0.98 / 0.99 | 3 / 9 / 46 / 125 | 0.933 / 0.964 / 0.989 / 0.993 | 0.929 / 0.962 / 0.986 / 0.993 | 0.979 / 0.976 / 0.990 / 0.993 | 0.973 / 0.969 / 0.983 / 0.990 |

Reading: 15 SubFam rows contain the AluY, AluYa5 and AluYb8 consensuses exactly and AluYb9 at 0.993; VSEARCH needs about 140 rows to reach a similar identity
(0.993-0.997) and CD-HIT-EST 125. A centroid is one copy and carries that copy's private mutations; a chunk consensus averages them away, which is the point of the method.
Limits: 316 copies, Y lineage only, one run per setting; the four references are Price's consensuses, the same ones used to label these copies, so this illustrates the
metric and the mechanism, it is not a validation. The independent test is the hg38 protocol in ../alu_hg38/. Tool versions: VSEARCH 2.27.0; CD-HIT-EST from Ubuntu's `cd-hit` package (parameters -n 8 at 0.90, -n 10 otherwise).
