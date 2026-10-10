# Real L1 copies (hg38 chr4, L1PA2-L1PA8): SubFam, peel into groups, group consensuses

Run 2026-10-08 in the cloud session (4 cores). Data: UCSC hg38 `rmsk.txt.gz` and `chr4.fa.gz` (hgdownload.soe.ucsc.edu). Not a benchmark with a truth; a worked example on real copies.

1. `extract_l1.py rmsk.txt.gz chr4.fa.gz l1.fa 1000`: every RepeatMasker L1PA2-L1PA8 copy on chr4 of at least 1 kb, without N, on the element strand. 2,222 copies (L1PA2 179, PA3 428, PA4 393, PA5 372, PA6 220, PA7 483, PA8 147); mean length about 3.3 kb.
2. `SubFam.sh -n 20 -c -t 4 -x l1 l1.fa`: 111 chunk consensuses, 5 min 2 s (`subfam.log`).
3. `peel_features.py` of the SINE-discriminator repository on the chunk consensuses (default settings, `peel.log`): 7 groups over 3 rounds, 105 of 111 chunks, 6 left over (G0).
4. `build_groups.py`: for every group, the consensus of its ORIGINAL copies (a random sample of at most 60, coverage-relative plurality 36 %, MAFFT `--retree 2`), not a consensus of chunk consensuses; then one MAFFT alignment of group heads and SubFam rows, grouped. Output: `examples/real_l1pa_chr4_n20c.aln.fasta`.

Row names: `CONS_G<k>_<sampled>of<all>copies` = group head; `<chunk>|G<k>|<rmsk majority>|p<purity>` = SubFam chunk consensus, the RepeatMasker subfamily that most of its copies carry, and the fraction of its copies carrying it.

## What it shows (read before quoting)
- 105 SubFam rows against the head of their own group: median identity 0.991 (min 0.812, max 1.000) over the shared span, median span 2.6 kb. For 30 of the 105 rows another group's head is marginally closer over the row's span: the heads are not far apart.
- The groups do not reproduce the RepeatMasker labels: chunk purity against rmsk is 0.69 on average; groups G1, G4, G5 mix 3-4 RepeatMasker subfamilies (e.g. G4: PA4 286, PA5 267, PA6 95, PA3 33 of 700 copies); G3 is mostly PA7/PA8, G2 and G7 mostly PA2/PA3. L1PA2-8 are age classes along one lineage, and the rmsk label is itself the best score to a consensus; there is no independent truth here.
- Peel group sizes are uneven (35 chunks vs 5) and G4 was peeled on only 4 features, so it may be several lineages.
- The sample of at most 60 copies per head means G4's head rests on 60 of 700 copies.
Neither "SubFam recovers the L1PA subfamilies" nor the opposite follows from this. It shows what SubFam's output and the peel give on real L1 copies, for the reader to inspect.


## 2. How the L1PA names were established, and why our groups disagree with them (research 2026-10-09)

**What was verified from sources (links), and what was not.**
- The classification is by the **3' end**. Smit, Toth, Riggs & Jurka 1995 (J Mol Biol 246:401-417, [PMID 7877164](https://pubmed.ncbi.nlm.nih.gov/7877164/)) analysed the 3' ends of about 900 primate L1 elements and found 47 subfamilies (full text read from the PDF supplied by the owner, pages cited below). Khan, Smit & Boissinot 2006 ([Genome Res 16:78, PMC1356131](https://pmc.ncbi.nlm.nih.gov/articles/PMC1356131/), full text read) say: over 90 % of L1 copies are 5' truncated, so the families were classified by their 3' ends, and copy numbers were counted as 3'UTRs per family.
- The families L1PA17 to L1PA1 are **one lineage**, each family replacing its predecessor (Khan 2006, Fig. 2; "a single lineage of L1 families amplified over the last 40 Myr"). L1PA2-8 are therefore consecutive time slices along one line of descent, not separate branches.
- Khan 2006 built full-length consensuses from RepeatMasker copies longer than 6 kb (April 2003 table), confirmed the family of each element by phylogeny and by comparing its 3'UTR with the Smit 1995 consensuses in Repbase, and removed CpG-decay mutations from the consensus. Ages (their Table 1): L1PA2 7.6 Myr, PA3 12.5, PA4 18.0, PA5 20.4, PA6 26.8, PA7 31.4, PA8 40.9; average divergence of copies from their family consensus: 1.24, 2.21, 3.02, 3.39, 4.44, 5.16, 6.92 %.
- The UCSC hg38 `rmsk` table (the one used here) was made in 2010 with RepeatMasker 3.0.1 and the Repbase library of 2010-03-02, and has not been recomputed ([UCSC track description](https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/hg/makeDb/trackDb/rmsk.html)). Repbase is not available to us; the **current** Dfam 4.0 holds L1PA2-L1PA8 only as **3'-end models of about 900 bp** (families `L1PA2_3end` ... `L1PA8_3end`, [DF000000339-345](https://www.dfam.org/family/DF000000340/summary), author list Finn, Hubley, Jones, Jurka, Smit, Wheeler; citation Smit 1995; "ORF2 ends at bp ~684").
- Not verified: the content of the 2010 Repbase library itself, how RepeatMasker 3.0.1 breaks ties between near-identical library members (code not read), and whether Dfam's 2026 consensuses equal the 2010 Repbase ones (they carry the Repbase name as alias, but the sequences may differ).

**Measured here (2026-10-09), with the Dfam 4.0 3'-end consensuses.**
1. The seven consensuses differ very little from their neighbours over the aligned ~900 bp (MAFFT, columns that differ incl. gaps): PA2/PA3 8, PA4/PA5 11, PA5/PA6 20, PA6/PA7 31, PA7/PA8 31; PA2/PA4 26, PA2/PA8 92 (identity 0.991, 0.989, 0.978, 0.966, 0.986 for the adjacent pairs). A single copy carries 1-5 % private divergence from its consensus (Khan 2006), i.e. 10-45 private differences in 900 bp, as many as or more than the 8-31 differences that separate adjacent families.
2. Re-labelling our 2,222 copies by the last 900 bp (best hit, `relabel_3end.sh`; 1,921 copies have a hit at 50 % identity or more, 301 have none, probably because they lack the 3' end): the label equals the hg38 rmsk label for 1,760 of 1,921 (91.6 %), is within one class for 1,843 (95.9 %), and every disagreement is between neighbouring classes (e.g. rmsk PA3 -> PA2 25, rmsk PA5 -> PA4 27, rmsk PA6 -> PA5 16). 72 copies tie between two consensuses. So the rmsk label is reproducible from the 3' end, and its remaining noise sits exactly at the neighbouring-class boundaries.
3. Our peel groups combine neighbouring classes (3'-end re-label of the copies in each group): G2 = PA3 123 + PA2 103; G7 = PA3 90 + PA2 13; G6 = PA2 52 + PA3 31; G3 = PA7 362 + PA8 99 + PA6 33; G4 = PA4 301 + PA5 260 + PA6 79; G5 = PA6 43 + PA5 22 + PA7 14; G1 (14 chunks, peeled on 76 features) is a mix of all classes. No group joins distant classes, such as PA2 with PA7.

**Reading.**
- The L1PA names are age classes along one lineage, defined by a few dozen diagnostic positions in the last ~900 bp, in a library built by Smit (1995) and refined by Jurka and Smit. A copy gets the name of its best-scoring consensus, which for a 5'-truncated copy depends on that 3' end only.
- Our groups come from the whole copy (several kb) and from what the peel finds shared between chunk consensuses; they merge neighbouring classes (PA2+3, PA4+5(+6), PA7+8) and never distant ones. The disagreement is one of granularity, and it coincides with the places where the consensuses themselves are only 8-31 positions apart, below the private divergence of a single copy.
- That is no evidence that the peel groups are wrong or that the rmsk classes are: both are partitions of a continuum. What would decide it is an independent test of whether PA2 and PA3 (8 differing positions) are separate lineages or one population; no source read here does that.
- For the paper: report that the L1PA labels are used only as an external age-ordered reference, not as ground truth, and state the 91.6 % reproducibility of the label from the 3' end as the noise floor of that reference.


### Smit et al. 1995, read in full: how a subfamily was defined (added 2026-10-09)
- **Data and method (p. 414-415, Materials and Methods).** The last 700 bp of the human L1Hs ORF2 were searched against the primate GenBank sequences (Release 79; 339 non-redundant ORF2 fragments longer than 50 bp, and 674 unique sequences with homology to the 3'UTR), aligned with CLUSTAL V and refined by hand, and a second route started from MER repeats. "To identify sequence subfamilies, we looked for multiple diagnostic differences between groups of homologous L1 sequences." A split was kept under one rule: **"Subfamily status was accepted when subdivision resulted in segregation of multiple diagnostic differences at separate positions in four or more sequences."** So a subfamily is a group of at least four database sequences that share several diagnostic differences, found by eye and by comparing alignments; there is no statistical test and no null model.
- **Consensus (p. 415).** The most abundant base per column; when two or more bases lie within 25 % of the top frequency, an IUPAC ambiguity code. Divergence excludes CpG sites. The authors note that "the continuous spectrum of evolving L1 source genes makes the division of L1 sequences into 47 subfamilies, by necessity, somewhat arbitrary" (p. 404), and that some consensuses contain ambiguous sites, "suggesting that the subfamily for which the consensus was derived may be split still further".
- **Order and names (p. 404-405).** The subfamilies are arranged by average divergence from their consensus (youngest first); the number counts from the subfamily derived from the most recently active source gene. The sample behind each subfamily is small: about 900 sequences for 47 subfamilies.
- **Two resolutions (p. 403, 407, Table 2).** The same sequences were also classified by ORF2 alone, giving ten subfamilies, and the authors report that the two classifications are "fully in accord, ... comprise sets of consecutive subfamilies". Table 2 gives the correspondence: ORF2-based L1P1 = **L1PA2-3**, L1P2 = **L1PA4-6**, L1P3 = **L1PA7-9**, L1P4 = L1PA10-16, L1P5 = L1PB1-3. The 3'UTR is the finer instrument "because the 3'-UTRs of L1 source genes vary to a much larger extent than the coding regions".
- **What this says about our groups.** The peel groups built from whole copies (which are mostly ORF-containing sequence, little 3'UTR) are G2 = PA2+PA3, G4 = PA4+PA5(+PA6), G3 = PA7+PA8(+PA6): almost exactly the coarser ORF2-level classes of Smit's Table 2 (PA2-3, PA4-6, PA7-9). The disagreement with the RepeatMasker names is therefore the difference between a classification by whole-copy sequence (ORF2 resolution) and one by the 3'UTR (the finer resolution that the names encode), as Smit already described in 1995. This is an observation of agreement in grouping, not a test that our groups are the ORF2 classes: the groups were not compared with the 1995 ORF2 consensuses.

## 3. Re-run with SubFam 1.3.0, grouped by the RepeatMasker label, and the 3'-end control (2026-10-10)

The peel-based grouping of section 1 is kept as a record but is no longer used for
`examples/real_l1pa_chr4_n20c.aln.fasta`. That file is now built like the other examples
(`examples/rebuild.sh`): `SubFam.sh -n 20 -c` on the same 2,222 copies (111 chunks; 1.3.0 gives
the same chunking as 1.2.0 here, 2,222 / 20 is a whole number of chunks either way), then for
each L1PA class the consensus of at most 60 randomly chosen member copies (seed 1) as the head row,
followed by the SubFam rows whose majority is that class.

**3'-end control.** The L1PA names are defined on the 3' end (Smit 1995, section 2). Running
SubFam `-n 20` on the last 900 bp of the same copies (`l1_3p.fa`, 48 s) and scoring chunk purity
against the rmsk labels:

| input | chunks | mean purity | within one class | chunks with purity >= 0.9 |
|---|---|---|---|---|
| full-length copies, `-c` | 111 | 0.696 | 0.926 | 32 |
| last 900 bp | 111 | 0.796 | 0.947 | 61 |

So the same copies sort by their rmsk class markedly better when only the region the classes were
defined on is used; the rest of the disagreement is of the order of the labels' own uncertainty
(re-labelling by the 3' end reproduces rmsk for 91.6 % of copies, section 2). Full-length copies
group at the coarser ORF2 resolution of Smit's Table 2; that is a property of the sequences, not an
error of either classification.
