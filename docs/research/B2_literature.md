# B2 SINE subfamily classification — literature and database report

Research note, 2026-10-10. Paper metadata and full texts via PubMed/PMC (DOIs given); database
facts from the Dfam REST API and the UCSC REST API, queried directly. Items marked **not verified**
were not opened. Compiled by a research pass; cross-check before citing.

## 1. Per-paper summary (all opened unless noted)

| Paper | Journal / DOI | Taxa | B2 subfamilies | Criteria / method | Data availability |
|---|---|---|---|---|---|
| Vassetzky, Borodulina, Ustyantsev, **Kosushkin, Kramerov** 2021 "Analysis of SINE Families B2, Dip, and Ves…" (full text read) | Int J Mol Sci 22:9897, https://doi.org/10.3390/ijms22189897 (PMC8466645) | Mouse mm10, rat Rnor_6.0, Chinese hamster, deer mouse, blind mole-rat (Spalax) | **5: B2a, B2b, B2c, B2d, B2e** (+ variants B2b1–3, B2c1–4 in cricetids; B2d1/d2 in Spalax) | 153,991 B2 hits in mm10 by custom Smith–Waterman Perl scripts (≥65% id, ≥90% length); 42,882 dimeric copies excluded; **92,352 full-length monomeric copies** split by "relatively long indels": B2a+B2b share an 11-bp central deletion; B2d has 28-bp and B2e 9-bp 3′ insertions; B2a differs in a 27-bp region before the tail (the τ signal). Subfamilies identified "manually and/or by a domestic script" on MAFFT alignments + fuzznuc/CD-HIT (id≈0.65); only subfamilies ≥1% of copies kept. Age by mean similarity of 100 random copies (B2a 83%, B2d 62%, B2e 59%) and TSD retention. Order (old→young): B2d, B2e → B2c, B2b → B2a. B2a = 64% of mouse B2; murid-only; absent in cricetids; Spalax has only B2d. Rat ~121,000 monomeric B2, same 5 subfamilies. | Consensuses only in Fig. 1A (figure). No GenBank/Repbase accessions, no alignment files mentioned in the text read (Supplementary section not retrievable; MDPI blocked). |
| Borodulina, Kosushkin, Ustyantsev, Vassetzky, Kramerov 2025 "Analysis of RNA Transcribed by Pol III from B2 SINEs" (full text read) | Noncoding RNA 11:39, https://doi.org/10.3390/ncrna11030039 (PMC12101331) | Mouse GRCm38.p6 | Same 5 (a–e); states B2a = "previously recognized as the four subfamilies Mm1a, Mm1t, Mm1o, and Mm2" | SSEARCH36 with B2a consensus (65%/90%), iterative depletion; copies assigned by best bitscore to subfamily consensus. ≥59,000 B2a in mouse, ≥85,000 in rat. B2 length "~180 nt". Transcript share: B2a 91%, b 4%, c 3%, d 2%, e 0.3%. | Supplementary worksheets list transcribed copies; consensuses not deposited (not stated). |
| Ichiyanagi T. … **Ichiyanagi K.** 2021 "B2 SINE copies serve as a transposable boundary…" (full text read) | Mol Biol Evol 38:2380, https://doi.org/10.1093/molbev/msab033 (PMC8136502) | Mus musculus B6 vs MSM | Repbase B2_Mm1a, B2_Mm1t, B2_Mm2 **+ new B2_Mm1o** ("older") | NJ tree (Clustal Omega/MEGA5) of 1,241 polymorphic + 6,000 random copies; Mm1o consensus built from an intermediate clade carrying both Mm1- and Mm2-specific nucleotides. Re-annotation with RepeatMasker custom library: **Mm1a ~17,000; Mm1t ~11,000; Mm1o ~22,000; Mm2 ~72,000** (total B2 ~120,000). Age: Mm2 oldest, Mm1o intermediate, Mm1a/Mm1t youngest. Northern: B2 RNA ~190 nt. | Mm1o consensus shown in Fig. 1C only; GEO for sequencing data. Not in Dfam (API search "Mm1o" = 0 hits). |
| Kawase & Ichiyanagi 2023 review "Mouse retrotransposons…" (HTML opened) | Genes Genet Syst 98, https://doi.org/10.1266/ggs.23-00221 | mm39 | Lists B2_Mm1a, Mm1t, Mm1o, Mm2, B3, B3A | From RepeatMasker mm39: B3/B3A divergence 23.3–26.8%; B2_Mm1a/1t/1o/Mm2 5.5–11.2%. B2 ≈190 bp, B3 ≈210 bp; first ~130 bp shared. Notes Mm1o absent from UCSC mm39 rmsk. | Supplementary Table S1 (averages). |
| Horton … Chuong 2023 "Mouse B2 SINE elements function as IFN-inducible enhancers" (full text read) | eLife 12:e82617, https://doi.org/10.7554/eLife.82617 | mm10, Dfam v2.0 annotation | B2_Mm2 (80,541), B2_Mm1a (16,321), B2_Mm1t (35,812) | Counts of annotated elements; consensuses from Repbase v24.02. | — |
| Kass, Kim, Rao, Deininger 1997 "Evolution of B2 repeats: the muroid explosion" (abstract only) | Genetica 99:1, https://doi.org/10.1007/BF02259494 | rat, mouse, deer mouse, hamster, gerbil, squirrel, hystricognaths | "distinct subfamilies" in mouse and rat; a "young mouse subfamily" matching transcribed B2 RNA | Copy numbers by hybridization: 80,000–100,000 per muroid genome; squirrel 2,500; <100 in guinea pig/nutria. 3′ end differs murid vs cricetid. N of subfamilies/diagnostic positions not in abstract — **not verified**. | not verified |
| Roy, Gong, Kass, Deininger 1998 (abstract) | DNA Seq 8:343, https://doi.org/10.3109/10425179809034078 | mouse | mentions "B2 type II consensus" | recent inserts match consensus | — |
| Bains & Temple-Smith 1989 (abstract) | J Mol Evol 28:191, https://doi.org/10.1007/BF02102476 | rodents | 71 B2 sequences; "seven regions of high variability… points to subfamilies" | by-eye/alignment of database sequences; no named subfamilies in abstract | **weak** (71 seqs) |
| Serdobova & Kramerov 1998 (abstract) | J Mol Evol 46, https://doi.org/10.1007/PL00006295 | Muridae, Cricetidae, Spalacidae; Dipodidae (DIP); Menetes (MEN) | No B2 subfamilies; defines B2 *superfamily* (B2, DIP, MEN, 4.5SI RNA; 80-bp core) | hybridization/PCR | — |
| Veniaminova, Vassetzky, Kramerov 2007 "B1 SINEs in different rodent families" (abstract) | Genomics 89:678, https://doi.org/10.1016/j.ygeno.2007.02.007 | 22 rodent families | **B1 only**, not B2 | indel-based variants | — |
| Kramerov & Vassetzky 2011 review "SINEs" (abstract) | WIREs RNA 2:772, https://doi.org/10.1002/wrna.91 | — | general review; B2 subfamily content **not verified** | — | — |

Note: PubMed finds no Kosushkin+Kramerov B2 paper before 2021. The 2021 IJMS paper is the
Kosushkin/Kramerov B2 subfamily paper.

## 2. Dfam / UCSC annotation (queried 2026-10-10)

| Dfam acc. | Name | Len | Clade | Citations in Dfam | Repbase alias | Consensus via API |
|---|---|---|---|---|---|---|
| DF000001735.1 | B2_Mm1a | 193 | Mus | Bains & Temple-Smith 1989; MGSC 2002 | B2_Mm1a | yes |
| DF000001736.1 | B2_Mm1t | 193 | Mus | same; description: "around 350000 genomic B2 copies in mouse" | B2_Mm1t | yes |
| DF000001737.1 | B2_Mm2 | 195 | Murinae | same | B2_Mm2 | yes |
| DF000001738.1 | B3 | 216 | Muroidea | Kalb et al. 1983 NAR | B3 | yes |
| DF000001739.1 | B3A | 198 | Muroidea | Kalb et al. 1983 | B3A | yes |

All five: classification `…SINE;tRNA_Promoter;No-core;L1-dependent;B2`, created 2015-09-08.
`https://www.dfam.org/api/families/<DF>/sequence?format=fasta` works. No rat- or hamster-specific
B2 entries; no B2_Mm1o.

UCSC rmsk, repFamily=B2, repName counts (mm39 whole genome; mm10 chr1 as cross-check):

| repName | mm39 count | median len (bp) | ≥170 bp | mm10 chr1 |
|---|---|---|---|---|
| B3 | 147,220 | 185 | 65% | 9,079 |
| B3A | 91,066 | 168 | 48% | 6,893 |
| B2_Mm2 | 87,172 | 176 | 59% | 5,647 |
| B2_Mm1t | 25,872 | 180 | 73% | 1,425 |
| B2_Mm1a | 18,416 | 188 | 81% | 1,035 |
| total | 369,746 | | | |

Consensus lengths 193–216 bp → full-length copies are ~180–200 bp.

## 3. Name mapping (Kramerov scheme ↔ Repbase/Dfam)

Stated explicitly in Vassetzky 2021 and Borodulina 2025:
- **B2a = B2_Mm1a + B2_Mm1t + B2_Mm2 (+ Ichiyanagi's B2_Mm1o)** — the "classical mouse B2", 64% of copies.
- **B2d = B3** ("originally described as B2l, B2 long"; the 28-bp insertion).
- B2b, B2c, B2e: no Repbase counterpart named; B3A not mapped in the text (**not verified**).
- The Kramerov group explicitly rejects Mm1/Mm2 as reliable subfamilies: using 13 diagnostic positions on
  5,000 B2a copies, only ~5% matched all Mm1 and ~14% all Mm2 characters; the rest are "hybrids"; they call
  Mm1/Mm2 "extreme cases in a continuum", though Mm2-like copies amplified earlier and Mm1-like later.

## 4. Solid vs weak

**Solid:** Dfam/Repbase has exactly 5 mouse B2-family entries with retrievable consensuses; B2_Mm1a/1t
differ by 2 nt and from Mm2 by 16 nt (Vassetzky 2021). Genome-wide counts above. The Kramerov
5-subfamily indel-based scheme is documented on 92k full-length mm10 copies with explicit diagnostic
indels and multi-genome support. Age ordering (B2d/B2e oldest → B2a youngest; within B2a, Mm2 → Mm1o →
Mm1) is consistent across the Ichiyanagi and Kramerov groups.

**Weak:** (i) Repbase Mm1a/Mm1t/Mm2 rest on single-nucleotide diagnostics with many intermediate copies;
Dfam cites only Bains 1989 (71 sequences) + MGSC 2002, no subfamily paper. (ii) B2_Mm1o exists only as a
figure consensus in one paper. (iii) Kramerov consensuses for B2a–e are not deposited anywhere found
(figure only); B2b/c/e have no Repbase name. (iv) Kass 1997 subfamily definitions unseen. (v) B3A's
position in the Kramerov scheme not stated.

## 5. Not verified / could not access
- Full text of Kass et al. 1997; Bains 1989 full text.
- Vassetzky 2021 Supplementary Materials (MDPI blocked); whether alignments/FASTA are supplied.
- Kramerov & Vassetzky reviews' B2 subfamily content (2005 Int Rev Cytol; 2011 WIREs).
- mm10 genome-wide counts (chr1 only).
- Any Ichiyanagi paper on B2 before 2021 (PubMed returned only the 2021 MBE paper).
- Repbase's own entry text/author lines.
