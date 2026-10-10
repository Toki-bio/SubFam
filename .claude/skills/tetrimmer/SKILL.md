---
name: tetrimmer
description: What TEtrimmer (TE library curation, Qian et al. 2025) does at the source level, its hard-coded limits (100 copies per MSA, 5 clusters, DBSCAN eps 0.1, genome required), and how to run or cite it in a comparison with SubFam. Use when a reviewer, a draft, or a benchmark mentions TEtrimmer, MCHelper, or "variant separation inside a TE family".
---

# TEtrimmer in one page

Read `docs/research/TETRIMMER.md` for the full profile with file:line references. The facts that
decide any comparison with SubFam, all verified in the source (commit 8954274, v1.7.4):

- Genome is **required** (`tetrimmer/TEtrimmer.py:100`); no entry point takes a FASTA of copies or
  an alignment.
- Copies per alignment: `--max_msa_lines 100`, chosen **longest first** (`functions.py:445-446`),
  random fill is unseeded.
- Clustering: divergent columns (major allele ≤ 0.8) → IQ-TREE K2P+I → patristic distance →
  `DBSCAN(eps=0.1, min_samples=3)` (`MSAcluster.py:366`), **not exposed** on the CLI; needs ≥ 20
  sequences, > max(50, 5 % L) divergent columns, clusters ≥ 10, largest ≥ 18 (`:584`), at most
  `--max_cluster_num 5`; noise usually dropped. Otherwise a single cluster with outliers removed.
- Consensus: A/C/G/T majority ≥ 0.7, gap-free. Cluster k of input X is named `X_0k`.
- Purpose: one boundary-correct consensus per library entry; the clustering protects the consensus
  from chimerism. It does not order or display the family and gives no assignment to most copies.

## How to word it in a paper

"TEtrimmer (Qian et al. 2025) splits the alignment of up to 100 copies of a library entry into at
most five DBSCAN clusters on a fixed tree-distance radius before building its consensus; the step
serves consensus curation and requires a genome. SubFam takes every copy of a family, uses no
distance radius or cluster count, and returns the family's variety as an ordered alignment."

## How to run it for a comparison

See `docs/research/TETRIMMER.md` §7: conda env from `TEtrimmer_env.yml` (needs Pfam-A, RepeatMasker,
IQ-TREE), seeds as `>name#SINE/B2`, `--debug` to keep the raw cluster BEDs, compare copy→cluster by
genomic overlap, align its consensuses into the SubFam example as `PUBLISHED_TEtrimmer_<k>` rows.
