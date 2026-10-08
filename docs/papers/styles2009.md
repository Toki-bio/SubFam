# Styles & Brookfield 2009 - source genes and gene conversion in the AluYh and AluYi lineages

BMC Evol Biol 9:102, doi 10.1186/1471-2148-9-102. Full text read (12 pp, PDF from the owner). Additional files 1-7 (alignments of all elements, publisher site) NOT read.

## What it is
A hand analysis of young Alu lineages. No algorithm and no significance test beyond a binomial estimate for parallel CpG mutation. Subfamilies = sets of copies sharing diagnostic mutations; "source genes" are inferred from groups of copies sharing extra mutations, strengthened when polymorphic (presence/absence) copies share them.

## Method (p. 10-11)
BLASTN of the human genome with AluYh9 and AluYi6 consensuses (1,426 and 1,277 hits), manual inspection for diagnostic mutations; not all diagnostic sites required because back mutation can restore ancestral bases. 1,000 bp flanks extracted, BLAT to chimpanzee (or human), orthologous regions aligned with ClustalW; polymorphism checked in the NCBI trace archives.

## Counts (Tables 1-3)
- AluYh7 (the 7-mutation core of AluYh9): 20 human-specific copies, 16 previously called AluYh9; Yh9 itself 2 copies; polymorphic.
- AluYh3a1 (4 diagnostic mutations): 98 human, 73 chimp, also gorilla/orangutan; AluYh3a3 (adds 19 bp deletion at 242-260): 3 human, 11 chimp.
- AluYi6: 237 human, 91 chimp; derivative putative subfamilies defined by 1-3 extra mutations (e.g. 57T 36 copies; 254A 35; 151T 77) -> up to 14 putative source lineages; several polymorphic.
- AluYg6: 380 human; Yg6a2 40, Yg5b3 27; at least 5 source elements.
- Yh7 and Yi6 appear to have multiple source genes; Yh3a3 is explained by a single master gene whose copy was inactivated by a backwards gene conversion in humans.

## Gene conversion and homoplasy (relevant to the benchmark)
- Complete conversion: a chimp AluYh3a1 element (DC7) sits at the orthologous position of a human AluSq; an old subfamily replaced by a young one (Fig. 3).
- Partial conversion: ancestral bases introduced at all four diagnostic sites of a 64 bp region (Fig. 4); orthologues share mutations outside the tract but not inside it. The authors say a conversion cannot be proven against parallel mutation in some cases.
- CpG transitions occur at about six times the rate of other mutations, causing homoplasy. Their calculation: in 98 Yh3a1 sequences the chance that at least one of 44 CpG sites mutates in >=14 sequences is >20 %.
- The 10-20 % figure in the introduction is the share of a young subfamily that "may operate as secondary source genes" (citing Cordaux et al. 2004), NOT a gene-conversion rate.
- Idea worth keeping: shared extra mutations among presence/absence-polymorphic copies cannot be old parallel decay, so polymorphism is an independent test of whether a diagnostic pattern is real. This addresses the homoplasy problem that the SINEderella manual parks.

## What it means for us
- Truth for these lineages is inferred and the paper states its own uncertainty (polymorphism data "cannot conclusively determine" absence; parallel mutation cannot always be excluded).
- Many derivative subfamilies (3-77 copies) are far below SubFam's default chunk of 50.
- Additional files give complete element alignments in human and chimp: http://www.biomedcentral.com/content/supplementary/1471-2148-9-102-S1.zip (and -S2 ... -S7), URLs as printed on p. 11.
