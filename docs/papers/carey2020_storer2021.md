# Notes: Carey et al. (Research Square preprint, 2020) and Storer et al. (Curr Protoc 2021)

Both read from PDFs supplied by the owner. Carey: all 17 pages. Storer: pages 1-20 of 40 (main protocol; the rest not read).

## Carey, Patterson, Wheeler - "TE subfamily annotation has a reproducibility problem"
NOT peer reviewed (preprint v3, 14 Dec 2020, CC BY). hg19 / hg38-panTro4, RepeatMasker open-4.0.5, Repbase 20140131.
- Replicates: TE pairs from single segmental duplications (16,962 filtered pairs) and human-chimp liftOver pairs. "Discordant" = RepeatMasker gives different subfamilies.
- Segmental duplications: Alu 10,347 pairs, 1,290 discordant (12.4 %); L1 6,615 pairs, 933 (14.1 %). Human-chimp: Alu 1,093,387 pairs 14.95 %; L1 1,050,856 pairs 17.60 %. By Alu type (Table 2): AluJ 16.8 %, AluS 12.2 %, AluY 15.5 % mismatch.
- Point-mutation drift model (diagnostic sites flipping) explains about half; ~7-8 % of AluS/AluY instances expected to look like another subfamily, mostly inside cliques (AluSg/Sx/Sx1/Sx3/Sz; AluYc/Yf1/Ym1/Y). AluJ unlikely to switch by point mutation.
- Apparent recombination (identity dips in 100 bp windows): discordant pairs show more (AluY 27.6 % vs 5.1 % concordant), but <30 % of discordant AluY pairs.
- Confidence from the ensemble of competing alignment scores (softmax of scores): median lower confidence 73.4 % in concordant vs 52.8 % in discordant pairs.
- Recommendations: use reliability measures to decide how aggressively to split families; annotation pipelines should report confidence.
- Relevance: the single-best-score "adjudication" criticised here is what a flat bitscore vote does. SINEderella's 10-cycle unanimity and the LEAK flag are a form of confidence; the paper gives a reference for the problem and a score-ensemble confidence measure to compare with.
- Data/scripts: http://wheelerlab.org/pubs/2020-discordant-CareyPatterson/CareyPatterson_suppl.tar.gz (blocked from the cloud container).

## Storer, Hubley, Rosen, Smit - Curation guidelines for de novo generated TE families
- Subfamily analysis uses COSEG: alignAndCallConsensus.pl -> bestwindow.pl -> preprocessAlignments.pl -> runcoseg.pl. Copies are aligned pairwise to ONE consensus (transitive MSA); stated to work when the reference "lack[s] significant indels".
- COSEG only uses sequences aligned over a chosen consensus window (bestwindow.pl). Example: consensus 246 bp, window 122 bp captures 98 of 117 alignments.
- Quote (p. 12): "In the event that the subfamily structure … presents with a 'truncated' appearance … COSEG is not appropriate because … a large section of sequences will be missed." For deletion/recombination products they suggest cd-hit-T + ClusterPartialMatchingSubs.pl (clusters by length then identity, 90 % then 80 %, hierarchically). It "works particularly well" for nearly identical copies (class II internal deletions, solo LTRs), "as opposed to an MSA with a wide range of divergences".
- COSEG options: -m minimum subfamily size (default 50), -u minimum distance between diagnostic sites (default 10; TG/CA pairs are one CpG decay, not co-segregation), -t three-site test.
- COSEG output can still be a mixture ("subfamily0 … no obvious pattern"); when COSEG finds nothing they group by divergence level. Visual assessment of the alignment (-html) is part of the protocol.
- Consensus calling is not a plain majority: highest-scoring nucleotide with a substitution matrix and CpG prediction.
- Relevance: Dfam's own curators state COSEG does not fit truncated/indel-structured sets and that visual review is part of the workflow.

## Where the Price 213 consensuses are (UNVERIFIED location)
Price et al. say the consensuses, full tree and subfamily list are in "Supplemental material … online at www.genome.org" (pp. 2245, 2247, 2249). Not reachable from the container. The COSEG repo's ALU.cons holds only AluSx.
