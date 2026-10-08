# Price, Eskin & Pevzner 2004 - notes from the full text

Genome Res 14:2245-2252, doi 10.1101/gr.2693004. Read in full (PDF supplied by the user,
pp. 2245-2252). The online supplement (consensuses of all 213 subfamilies, full MST) and the
authors' code (http://www.cs.ucsd.edu/~aprice/alu.html) have NOT been read.
The PDF itself is not committed (journal copyright).

## Data
- BLAST (bl2seq) of human genome Build 34 against the Repbase AluSx consensus (or equivalently RepeatMasker). 
- Multiple alignment = tabulate the nucleotide of each element at each position 1-282 of the AluSx
  consensus; insertions recorded separately.
- Elements missing >5 bases at either end are dropped: ~480,000 full-length Alus (of >1 M).
  Truncated copies are therefore NOT classified at all.

## Method (de novo, not seeded with Repbase)
Starts from ONE subfamily holding all ~480k elements and splits recursively. Repbase is only
used afterwards for comparison.
1. Why not k-means/EM (Table 1B): a nested small subfamily (Ya5, 0.7 %) is not found; k=2 splits 90/10 on a
   random mutation (pos. 2 G/T) instead. Single non-consensus frequencies are uninformative
   because per-position mutation rates are unknown and vary widely.
2. Biprofiles: look for over-represented PAIRS of non-consensus values. Ya5 example, 5 diagnostic
   positions (AluSx coordinates 91, 98, 146, 175, 238): C@1 with A@2 occurs 12x more often than
   expected from the marginals (0.024*0.023); P < 1e-2000. Their split gives 99.3 % / 0.7 %,
   the 0.7 % cluster is exactly Ya5 (C,A,A,T,C).
3. Linkage statistic for mutations mu1, mu2 in a subfamily of N elements: N12 observed vs N1*N2/N.
   P-value xi = hypergeometric-type tail: sum over N^ = N12..min(N1,N2) of
   multinomial(N; N^, N1-N^, N2-N^, N-N1-N2+N^) / (C(N,N1) * C(N,N2)).
   Subfamily xi = min over pairs, Bonferroni (x number of pairs tested). Threshold 0.001.
4. Safeguards (all of them shape what can be found):
   - Bonferroni checked on simulated uniform data: all P > 1.
   - Age bias: copies with mu1 tend to be older and also carry mu2; inflates N12 by up to 4/3; xi corrected.
   - Indels excluded; no two mutations closer than 10 nt.
   - CpG->TpG / CpG->CpA (and reverse) mutations excluded (non-independent, hypermutable).
   - Mutations to a value already present in other subfamilies excluded: to avoid assigning
     mosaic Alus (Alu-Alu recombination / gene conversion) to a new subfamily.
   - Pairs with N12 < 50 excluded -> minimum subfamily size 50.
   - After each split, all elements reassigned to the subfamily whose consensus fits, consensuses
     recomputed, iterate to convergence (assignment rule itself is not spelled out in the text).
   - Validation: the union of two or more subfamilies must FAIL the uniformity test, else merge
     (only checked for unions formed while building the MST, not for every possible union).
5. Result of pass 1: scaffold of 60 subfamilies.
6. Calibration: mutation rate at each position, consensus value -> each other value, averaged over
   elements of all scaffold subfamilies with that consensus value, relative to the subfamily's
   divergence from its consensus (proxy for age).
7. Pass 2: split each subfamily on over-represented INDIVIDUAL mutations (binomial test, P < 0.001,
   Bonferroni, same exclusions, min size 50; mutations present in an adjacent MST subfamily excluded
   unless that would just split off an intermediate). +153 subfamilies = 213 total.
   Sizes ~50 to ~60,000, most >= a few hundred. P-values <1e-6000 down to near 0.001, most < 1e-40.

## Tree
Subfamilies can be internal nodes (fossil record), so ordinary phylogeny does not apply.
Tree = Minimum Spanning Tree (Kruskal) on Hamming distance between consensuses (CpG transitions
ignored, indels penalised), rooted at the subfamily with the highest average divergence from its
own consensus. Number of leaves = lower bound on number of source elements: 143 for the 213
subfamilies (vs 14 for the 31 Repbase subfamilies). Age scale: 4 My per 1 % divergence
(AluJ ~15 % / 60 My, AluS 11-6 %, AluY 1-4 % / 4-16 My).

## Numbers worth keeping
- Repbase had 31 Alu subfamilies (+3 monomeric ancestors ~140 bp, outside the study).
- Their 213 include only 19 of the 31 Repbase subfamilies; the other 12 are minor AluY branches.
  213 + 12 = 225 subfamilies "previously and presently identified".
- Table 2 (size / P): AluJo 7,266 / 5e-11262 (corrected from the supplement; an earlier version of this line gave 8e-1841, which belongs to AluSx_3); AluSx 39,724 / 6e-4770; AluSq 4,035 / 2e-62;
  AluSp 28,063 / 7e-4520; AluY 27,023 / 2e-6924; AluYa5 3,257 / 4e-2813; AluSx_3 3,292 / 8e-1841;
  AluSx_5 401 / 3e-150; AluSq_3 1,956 / 2e-779; AluSg_4 1,904 / 1e-679; AluSc_8 9,588 / 1e-5959;
  AluY_8 107 / 1e-48. Repbase AluSq is smaller in their allocation (elements moved to novel neighbours).
- Their consensuses differ from Repbase mostly at CpG positions (ill-determined).
- Novel names (AluSx_3 etc.) are THEIR labels, given by attaching to the nearest Repbase name.
  They are not known to be the same as Dfam's AluSx1/AluSx3 (UNVERIFIED, do not equate).

## Supplement (SUBFAMILIES.FINAL.txt, supplied by the owner; checked 2026-10-08)
Not committed (publisher file; kept outside the repo). Verified by parsing:
- 213 subfamilies, ids 0-212, no duplicate names; 60 without a parent P-value (the scaffold), 153 with one (pass 2). Sizes sum to 476,152 copies; smallest 53, largest 58,050 (AluSz); consistent with the 50-copy minimum.
- Consensus lengths 279-288 (mostly 280-283); consensuses are therefore not all one length (they carry indels relative to AluSx).
- Matches the paper's Table 2 for AluSx 39,724; AluY 27,023; AluYa5 3,257; AluYb8 2,285; AluSx_3 3,292; AluSc_8 9,588; AluY_8 107; AluSx_5 401.
- The 12 Repbase subfamilies not found: AluYa8, Yc2, Yd2, Yd3, Yd3a1, Yd8, Ya4, Yf1, Yh9, Ybc3a, Ye2, Yf2.
- The six novel names with underscores shown in Table 2: Sx_3 (id 12), Sq_3 (15), Sx_5 (21), Y_8 (30), Sg_4 (53), Sc_8 (56).
- Tree image supplied (appears to be the minimum spanning tree over these ids; roots of the three boxes: 60 AluJo, 0 AluSx, 1 AluY). The caption/legend was not supplied, so the colour meaning is inferred: pink nodes are exactly the six novel subfamilies above.

## Authors' stated limitations
Excludes indels, frequent CpG mutations, and mutations already present in other subfamilies, so
subfamilies defined by those are hard to find. The partition need not be unique; subfamilies
whose elements are spread over more than one member of the partition may exist. No fix offered.
They conjecture thousands of source elements and many unidentified subfamilies; most Alus are
not retrotransposition-competent.

## What this means for SubFam
- Closest published goal to SubFam for Alu: a full-length, threshold-free split of ~480k copies.
- Their method is fully specified in the Methods and can be re-implemented as an independent,
  statistically defined reference partition on the same copies (needs a collinear alignment to a
  single reference: Alu qualifies, truncated copies are dropped).
- It is a defensible reference partition, not ground truth: partition not unique, parallel
  mutation and mosaic copies are deliberately excluded.
- Not stated in the paper: how well individual copies are assigned (no error rate).
