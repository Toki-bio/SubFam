# Ahmed, Li, Liang 2013 - three new Alu Yb subfamilies by source tracking

Mobile DNA 4:25, doi 10.1186/1759-8753-4-25. Full text read (11 pp, PDF from the owner). Additional files (Table S1 = insertions in the reference genome, Table S2 = outside it, Figures S1-S4) NOT read.

## Not a statistical method
Subfamilies are defined by diagnostic mutations found by inspection of an alignment, then ordered by nesting. No significance test; Price 2004 is cited only (ref. 33) as one of the sources of doubt about the master-gene model, not used.

## Data and method
- 3,179 Yb copies in the human genome, 80 % (2,545) human-specific (hsYb; list from an unpublished study, Tang et al.). Only full-length hsYb copies used; 1,322 hsYb8 studied.
- Source tracking: each hsYb copy BLASTed (e-value 1e-5) to the reference genome; copies whose best match was not unique were dropped (source cannot be determined); the rest aligned to Yb7/8/9 consensuses, grouped into similarity clusters, one neighbour-joining tree per cluster rooted on the Yb8 consensus; median-joining networks for some. Poly-A removed.
- Results: one Yb8 copy at hg19 chr10:10493416-10493732 heads a cluster of 714 Yb8 elements, said to have produced 713 copies (abstract: 712) = 54 % of 1,322 hsYb8; eight other Yb8 copies produced >=10 copies each ("stealth drivers"). One Yb8 copy at chr14:101990881-101991202 is the likely source of all hsYb9 (131 full-length Yb9 copies, 16 Yb8 best matches).
- New subfamilies, by diagnostic changes relative to AluY numbering: Yb8a1 = 259 G>A; Yb10 = 259 G>A + 174 C>G (the Yb9 mutation); Yb11 = 259 G>A + 174 C>G + 200 +T insertion. Found by noticing that every copy with the 200 +T insertion also carries the other two, whereas the reverse does not hold, hence a chronological nesting.
- Copies found: reference genome Yb8a1 75, Yb10 8, Yb11 16; 1000 Genomes / HuRef extra 6, 3, 15; totals 81, 11, 31. Polymorphic: Yb11 28/31 (~90 %), Yb10 5/11 (45 %), Yb8a1 28/81 (34 %).
- Validation of Yb11: 130 trace reads carry the T insertion; 3 of 15 outside-reference loci confirmed in trace data, 5 loci amplified by PCR and sequenced.
- Ages from non-CpG mutation density with a neutral rate 0.15 %/Myr: Yb8a1 1.93 Myr, Yb9 2.15, Yb10 1.73, Yb11 0.71. NOTE internal inconsistency: the abstract says Yb10 1.4 Myr, the text 1.73 Myr.
- Pathway proposed: Yb8 -> Yb8a1 -> Yb10 -> Yb11 (linear). Alternative for Yb10 (Yb9 gaining 259 G>A, perhaps by gene conversion) is "less convincing" but "requires further analysis". Molecular-clock test rejects uniform rate for Yb8a1 and Yb10, not for Yb11 (P < 0.43).

## What it means for us
- A small, nested, young truth set with explicit counts and single-nucleotide diagnostics: Yb11 (16 copies) and Yb10 (8) are far below SubFam's default chunk of 50 (would need -n well under 20), and Price p. 2246 states the pairwise approach "cannot identify subfamilies that differ at only a single diagnostic position".
- The labels are the authors' own diagnostic definitions, not independently derived; only Yb11 has direct sequence validation. Ambiguity about Yb10's origin is stated by the authors.
- Copy lists exist in Additional file 1 / 3 (publisher site, blocked here); could serve as a benchmark truth set for the Yb lineage once obtained. Names Yb8a1/Yb10/Yb11 not checked against Dfam (UNVERIFIED).
