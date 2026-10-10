# Related work: tools that overlap with SubFam, and comparison feasibility

Research note, 2026-10-10. "Unverified" = detail not confirmed against the paper/docs; DOIs marked
(mem) are from memory. Cross-check before citing.

SubFam: k-mer UPGMA order → fixed-N chunks → MAFFT + plurality consensus → short "proxy" alignment;
no identity threshold, no subfamily claim.

## A. TE subfamily inference

**COSEG** (Smit/Hubley/Price; RepeatMasker package). Finds subfamilies by statistically significant
co-segregating 2–3 bp mutations in near-full-length copies aligned to one consensus (cross_match/RMBlast
input, `preprocessAlignments.pl`); outputs subfamily consensuses and a GraphViz tree; README says its
consensus builder is "rather crude", optional `refineConsSeqs.pl` uses RepeatModeler's Refiner.
Shares: input = many copies of one family, output = several consensuses. Differs: explicit subfamily
claim via diagnostic positions + p-values; needs one reference consensus and full-length copies.
Availability: https://github.com/rmhubley/coseg and https://repeatmasker.org/COSEGDownload.html,
v0.2.4 (2024-07-30). Comparison: **medium** (build cross_match alignments of copies to a seed; `-m 50`
min-subfamily size is the analogue of `-n 50`). Measure: how many COSEG subfamily consensuses are
recovered (≥X% identity) by a SubFam row and vice versa; runtime; behaviour on truncated copies.

**Price, Eskin & Pevzner 2004** (Genome Res 14:2245, doi:10.1101/gr.2693004). Original Alu-subfamily
algorithm (213 Alu subfamilies vs 31 in Repbase) — COSEG descends from it, so citing both and
comparing with COSEG covers it. Original code URL likely dead (unverified).

**AnTE / AnTE2** (Wacholder … Pollock, PLoS Genet 2014, doi:10.1371/journal.pgen.1004482; AnTE2 bioRxiv
2016, doi:10.1101/094326). Bayesian ancestral-source inference; argues COSEG subfamilies lump multiple
source elements. Comparison: **hard** (code availability unclear); cite as conceptual contrast.

**RepeatModeler2 / Refiner** (Flynn et al. 2020 PNAS 117:9451, doi:10.1073/pnas.1921046117). De novo
family discovery; consensus by majority rule from a multiple alignment of sampled copies; families within
20% similarity labelled "putative subfamily". One consensus per family. Comparison: **medium**
("component" comparison). https://www.repeatmasker.org/RepeatModeler (OSL-2.1).

**MCHelper** (Orozco-Arias et al. 2024 Genome Res 34:2256, doi:10.1101/gr.278821.123) and **TEtrimmer**
(Qian et al. 2025 Nat Commun, doi:10.1038/s41467-025-63889-y; https://github.com/qjiangzhao/TEtrimmer).
Automated library curation; TEtrimmer separates variants inside one MSA via tree + DBSCAN clustering
before re-building consensuses — the closest recent idea to SubFam. Differs: requires a genome, starts
from one consensus, clustering uses distance parameters (DBSCAN eps). Comparison: **medium**.
**TE-Aid** (Goubert et al. 2022 Mobile DNA, doi:10.1186/s13100-021-00259-7): diagnostic plots, cite only.

**Classifiers, not overlapping** (cite to pre-empt): TEsorter (doi:10.1093/hr/uhac017), DeepTE (Yan et
al. 2020 Bioinformatics, doi (mem):10.1093/bioinformatics/btaa519), TEclass2 (doi:10.3390/biology15010059),
REPCLASS (unverified). **TE-greedy-nester** (Lexa et al. 2020, doi:10.1093/bioinformatics/btaa632):
nested LTR copies — unrelated. **pantera** (Sierra & Durbin 2024 Mobile DNA,
doi:10.1186/s13100-024-00323-y): clusters pangenome insertion polymorphisms at high identity — threshold-
based, cite as contrast. No dedicated 2023–2026 TE-subfamily-inference tool found beyond these
(unverified absence).

## B. Consensus from many copies

RepeatScout (Price et al. 2005 Bioinformatics 21:i351), EDTA (Ou et al. 2019 Genome Biol), Dfam seed
alignments → consensus/HMM, EMBOSS `cons` (SubFam reimplements it byte-identically), hmmemit. All give one
consensus per family; none order-and-chunk. Comparison: **easy** component test — single family-wide
consensus vs SubFam's k rows: per-copy identity to the nearest SubFam row vs to the single consensus.
MAFFT large-N modes (`--parttree`, `--retree 0 --reorder`, used by `-m`) are aligners, not reducers.

## C. Threshold clustering and ASV tools

CD-HIT (Fu et al. 2012, doi:10.1093/bioinformatics/bts565), VSEARCH (Rognes et al. 2016 PeerJ,
doi:10.7717/peerj.2584), MMseqs2/Linclust (Steinegger & Söding 2018 Nat Commun,
doi:10.1038/s41467-018-04964-5). Active, permissive licenses, trivial to run (`--cluster_fast --id X
--consout`). Comparison: **easy** and the one every reviewer will ask for. Measure: number of
clusters/singletons vs threshold across young vs old families, and how many ground-truth subfamily
consensuses are recovered.
**Swarm** (Mahé et al. 2022 Bioinformatics 38:267, doi:10.1093/bioinformatics/btab493) — threshold-free
local linkage d=1 with abundance-based breaking. Designed for amplicons where identical reads are
abundant; TE copies are nearly all unique, so d=1 chains collapse or fragment (**easy** short negative
demonstration). **DADA2** (doi:10.1038/nmeth.3869), **Deblur** (doi:10.1128/mSystems.00191-16): error-
model-based ASV inference; not applicable (divergence is mutation, not sequencing error) — cite only.

## D. Subset selection and divide-and-conquer MSA

**PASTA** (Mirarab et al. 2015 J Comput Biol), **UPP** (Nguyen et al. 2015 Genome Biol 16:124,
doi:10.1186/s13059-015-0688-z), **MAGUS** (Smirnov & Warnow 2021 Bioinformatics 37:1666,
doi:10.1093/bioinformatics/btaa992), **Recursive MAGUS** (Smirnov 2021 PLoS Comput Biol,
doi:10.1371/journal.pcbi.1008950). Divide via a guide tree into disjoint subsets, align subsets (MAFFT),
then *merge into one full-size alignment*. Key difference a reviewer will probe: SubFam collapses each
subset to a consensus and never merges copies — output is k rows, not n rows; the goal is reduction, not
full-MSA accuracy. Comparison: **medium**; a paragraph (merge vs collapse) plus a runtime/row-count line
suffices.
**Hobohm** (Hobohm et al. 1992 Protein Sci; https://github.com/agormp/hobohm), **Treemmer** (Menardo et
al. 2018 BMC Bioinformatics), **PARNAS** (bioRxiv 2022, doi:10.1101/2022.09.12.507613). Pick *existing*
representative sequences (need a cutoff or a tree); SubFam emits *averaged* sequences that remove private
mutations. Comparison: **easy-medium**; measure identity to ground-truth source sequences (consensus
should win for old families).
**MAFFT --add/--addfragments/--keeplength** (Katoh & Frith 2012, doi:10.1093/bioinformatics/bts578):
places copies back onto the proxy alignment; complementary, not competing.

## E. Ordering for visual overview

Jalview "Sort by tree" (https://jalview.org/help/html/calculations/sorting.html), MAFFT `--reorder`, R
`seriation` (Hahsler 2008 JSS). No published "k-mer UPGMA with leaf-flipping" ordering tool found besides
ViewAlign (SubFam's own source); cite seriation/optimal-leaf-ordering (Bar-Joseph 2001, unverified) as
prior art. Comparison: **not needed**; a note that ordering was checked against MAFFT guide-tree order
(`-m`) pre-empts it.

## F. Satellite monomer tools

**TAREAN** (Novák et al. 2017 NAR 45:e111, doi:10.1093/nar/gkx257): k-mer-graph consensus from reads.
**StringDecomposer** (Dvorkina et al. 2020, doi:10.1093/bioinformatics/btaa454), **CentromereArchitect**
(2021, doi:10.1093/bioinformatics/btab265), **HORmon** (Kunyavskaya et al. 2022 Genome Res 32:1137):
monomer inference + HOR annotation from an assembled centromere; need a template monomer.
**HumAS-HMMER** (Uralsky et al. 2019 Data in Brief; PMC6447721) and **Alpha-CENTAURI** (Sevim et al.
2016; PMC4920115): classify monomers against known classes. Comparison: **medium** for HORmon's monomer-
inference module (only if the paper keeps a satellite use case).

## Ranking: comparisons a reviewer will most likely demand

1. **VSEARCH/CD-HIT threshold sweep** (easy, essential).
2. **COSEG** (medium, essential for TE readers).
3. **PASTA/MAGUS** (conceptual rebuttal: merge vs collapse).
4. **Single family consensus** (easy): what one consensus hides.
5. **Representative subsampling (Treemmer/PARNAS/Hobohm)** (easy-medium).
6. **TEtrimmer cluster-and-consensus stage** (medium): nearest recent idea.
7. **HORmon** (medium): only with a satellite example.
8. Swarm (easy, optional negative demonstration).

## Minimal baseline set (argued)

Keep: (1) VSEARCH at 3–4 thresholds, (2) COSEG, (3) one family-wide consensus. PASTA/MAGUS: a paragraph.
Drop ASV tools, classifiers, TE-greedy-nester, TAREAN. TEtrimmer: mention; running its clustering step is
cheap insurance.

## References (verified unless marked)

- COSEG: https://github.com/rmhubley/coseg ; https://repeatmasker.org/COSEGDownload.html (v0.2.4, 2024-07-30)
- Price, Eskin, Pevzner 2004 Genome Res 14:2245. doi:10.1101/gr.2693004
- Wacholder et al. 2014 PLoS Genet 10:e1004482. doi:10.1371/journal.pgen.1004482 ; AnTE2 bioRxiv doi:10.1101/094326
- Flynn et al. 2020 PNAS 117:9451. doi:10.1073/pnas.1921046117
- Orozco-Arias et al. 2024 Genome Res 34:2256. doi:10.1101/gr.278821.123
- Qian et al. 2025 Nat Commun. doi:10.1038/s41467-025-63889-y ; github.com/qjiangzhao/TEtrimmer
- Goubert et al. 2022 Mobile DNA. doi:10.1186/s13100-021-00259-7
- Sierra & Durbin 2024 Mobile DNA 15:13. doi:10.1186/s13100-024-00323-y
- Lexa et al. 2020 Bioinformatics. doi:10.1093/bioinformatics/btaa632
- Zhang et al. 2022 Hortic Res. doi:10.1093/hr/uhac017 (TEsorter); Bickmann et al. 2025 Biology. doi:10.3390/biology15010059 (TEclass2); DeepTE doi (mem):10.1093/bioinformatics/btaa519
- Price et al. 2005 Bioinformatics 21:i351 (RepeatScout; DOI unverified); Ou et al. 2019 Genome Biol (EDTA; DOI unverified)
- Fu et al. 2012 Bioinformatics 28:3150. doi:10.1093/bioinformatics/bts565
- Rognes et al. 2016 PeerJ 4:e2584. doi:10.7717/peerj.2584
- Steinegger & Söding 2018 Nat Commun 9:2542. doi:10.1038/s41467-018-04964-5
- Mahé et al. 2022 Bioinformatics 38:267. doi:10.1093/bioinformatics/btab493
- Callahan et al. 2016 Nat Methods 13:581. doi:10.1038/nmeth.3869 ; Amir et al. 2017 mSystems 2:e00191-16. doi:10.1128/mSystems.00191-16
- Mirarab et al. 2015 J Comput Biol (PASTA; venue unverified); Nguyen et al. 2015 Genome Biol 16:124. doi:10.1186/s13059-015-0688-z
- Smirnov & Warnow 2021 Bioinformatics 37:1666. doi:10.1093/bioinformatics/btaa992 ; Smirnov 2021 PLoS Comput Biol 17:e1008950. doi:10.1371/journal.pcbi.1008950
- Hobohm et al. 1992 Protein Sci 1:409 ; github.com/agormp/hobohm
- Menardo et al. 2018 BMC Bioinformatics (Treemmer); PARNAS bioRxiv doi:10.1101/2022.09.12.507613
- Katoh & Frith 2012 Bioinformatics 28:3144. doi:10.1093/bioinformatics/bts578
- Hahsler et al. 2008 J Stat Softw 25(3) (seriation); Jalview sorting docs: https://jalview.org/help/html/calculations/sorting.html
- Novák et al. 2017 NAR 45:e111. doi:10.1093/nar/gkx257
- Dvorkina et al. 2020 Bioinformatics 36:i93. doi:10.1093/bioinformatics/btaa454 ; Dvorkina et al. 2021 Bioinformatics 37:i196. doi:10.1093/bioinformatics/btab265 ; HORmon Genome Res 2022 32:1137
- Uralsky et al. 2019 Data in Brief (HumAS-HMMER; PMC6447721); Alpha-CENTAURI: PMC4920115
