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
