# Scaling SubFam: beyond 30,000 copies, and to LINEs

Measurements on 4 cores, simulated families (`benchmark/simulate.py`, 300 bp copies unless
noted). "Purity" is the fraction of copies whose chunk majority is their own subfamily, computed
from the ordering alone (`-n 50`), so ordering methods can be compared without aligning anything.

## 1. Where the 30,000 wall is

SINEderella's copy of SubFam orders the copies with `mafft --retree 0 --reorder`. Without
`--parttree`, MAFFT computes the full all-against-all 6-mer distance matrix: O(n²) time and
memory.

| copies | `mafft --retree 0` | `mafft --retree 0 --parttree` | SubFam k-mer tree (`-k 6`, numpy) |
|---|---|---|---|
| 2,000 | 1.1 s, 22 MB | – | 1.7 s, 144 MB |
| 10,000 | 42 s, 403 MB | – | 27 s, 795 MB |
| 20,000 | 219 s, 1.6 GB | 60 s, 203 MB | 105 s, 2.2 GB |
| 100,000 | (extrapolated: ~1.5 h, ~40 GB) | 376 s, 1.1 GB | (n² matrix: 40 GB) |

So 30,000 is simply where the n² matrix stops fitting a workstation. SubFam 1.1 already has the
fix: **`-P`** (MAFFT PartTree, roughly O(n log n)). At 20,000 copies the PartTree order is as
good as the full one (8/8 subfamilies, purity 0.918 vs 0.930 for the k-mer tree). The SubFam
k-mer tree has the same n² limit as `--retree 0` and should be used below ~20,000 copies, where
it is the best ordering (purity 0.71 vs 0.60 on the old family).

Whole pipeline with `-P -n 50`:

| copies | ordering | chunk alignments + consensus | final L-INS-i of the consensuses | total |
|---|---|---|---|---|
| 20,000 (400 consensuses) | 61 s | 112 s | 47 s | 3 m 40 s |
| 100,000 (2,000 consensuses) | 6 m 14 s | 6 m 12 s | 17 m 35 s | 30 m 4 s |

Chunk alignment is linear and embarrassingly parallel. The ordering is no longer the bottleneck.
The 100,000-copy run recovered all 8 subfamilies with purity 0.944 and consensuses at 100 %
identity to the true sources.

### 1.1 The next wall: the final alignment

L-INS-i is O(m²·L²) in the number m of consensuses. 400 consensuses take 47 s, 2,000 take
17.5 min (22× for 5× the rows); 20,000 (a million copies) would take days. Two ways out, both cheap:

- **Hierarchical SubFam.** Treat the consensuses as a new input and run SubFam again
  (`-n 20`–`50`). Consensuses of consensuses converge to the subfamilies; `chunks.tsv` of both
  levels give every copy's lineage. This is what a human does with 2,000 rows anyway.
- **Align per partition.** When the input is partitioned first (§3), each partition's
  consensuses are a few dozen to a few hundred rows: L-INS-i stays fast, and the
  per-partition alignments are what one browses.

## 2. What VSEARCH does (code level) and what to borrow

Read from the VSEARCH source (`src/core/dbindex.cpp`, `unique.cpp`, `searchcore.cpp`,
`cluster.cpp`, `align_simd.cpp`, `msa.cpp`). Its clustering is **greedy single-pass centroid
clustering**, and everything is organised so that the per-query cost is linear in the number of
*centroids*, never in the number of input sequences:

1. **Index.** Each sequence contributes its set of *distinct* 8-mers (2-bit packed; k-mers
   overlapping non-ACGT or dust-masked bases are dropped). A k-mer present in ≥ 1/8 of the
   sequences gets a bitmap (one bit per sequence); rarer ones get a posting list. During
   clustering a sequence enters the index only when it becomes a new centroid.
2. **Candidates.** For a query, one 16-bit counter per indexed sequence is incremented for
   every shared k-mer (SIMD over bitmaps, scalar over posting lists, in cache-sized slices). Every
   indexed sequence with ≥ `minwordmatches` (12 for k = 8) shared k-mers is offered to a bounded
   min-heap of size `maxaccepts + maxrejects + 8`. Candidates are then taken best-first by
   shared k-mer count. `--id` plays no part in candidate selection.
3. **Verification.** Up to `maxaccepts + maxrejects` candidates (1 + 8 for `--cluster_fast`,
   1 + 32 otherwise) are aligned with a striped 16-bit SIMD global aligner, eight at a time
   (match +2, mismatch −4, gap open 20 / extend 2, terminal gaps 2 / 1). Identity is
   `--iddef 2` by default (terminal gaps excluded). The first candidate reaching `--id` is
   accepted (`maxaccepts 1`: the first acceptable centroid, not the closest); after
   `maxrejects` failures the query becomes a new centroid and is added to the index.
4. **Order.** `--cluster_fast` sorts by length (longest first), `--cluster_size` by abundance,
   `--cluster_smallmem` keeps input order. Order matters because centroids never move and later
   sequences only see earlier centroids.
5. **`--consout`.** No multiple aligner: each member's CIGAR against its centroid is stacked into
   a star alignment; per column the most frequent base wins if its count ≥ the gap count,
   otherwise the column is dropped.
6. **Long sequences.** The SIMD aligner handles pairs up to 25 M cells (5 kb × 5 kb); larger
   pairs fall back to a scalar Hirschberg aligner. A 6 kb LINE copy against a 6 kb centroid
   always takes the slow path. Shared k-mer counts are absolute (not length-normalised), so long
   sequences crowd short centroids out of the top-N.

**Why it cannot replace SubFam.** Candidate ranking is by exact shared 8-mers; a 8-mer survives
with probability ≈ identity⁸ (10 % at 75 %, 6 % at 70 %). Below ~75 % identity the true
neighbour often ranks outside the top `maxrejects` and is never aligned, and `--id` has to be
chosen per family (README benchmark: the working `--id` is 0.90 for a young family, 0.80 for
a middle one and exactly 0.70 for an old one).

**What to borrow.**

- *The centroid trick* for partitioning: a loose `--cluster_fast` pass is linear-ish and
  splits a million copies into groups SubFam can handle one at a time (§3).
- *The counter scan* as a nearest-neighbour finder: counting shared k-mers against an index is
  how to find each copy's top-m neighbours without a full matrix. A sparse m-nearest-neighbour
  graph plus single/average linkage on its edges gives an approximate guide tree in
  O(n·m·log n) memory, replacing the n² matrix of the k-mer UPGMA. This is the "fast k-mer
  approach" for SubFam proper; not implemented yet.
- *Iddef 2 and terminal-gap trimming* for comparing truncated copies (§4).

## 3. Partition first, then SubFam: measured

`vsearch --cluster_fast --id X` on the whole input, clusters with ≥ N copies become partitions,
all smaller clusters are pooled into one "rest" partition, SubFam (`-n 50`, k-mer tree) runs
inside every partition, and the consensuses are pooled (`benchmark/partition_test.sh`):

| family | `--id` | partitions (sizes) | subfamilies recovered | purity |
|---|---|---|---|---|
| middle (8 %) | none (plain SubFam) | 1 | 8/8 | 0.908 |
| middle | 0.80 | rest:1491 + 5 small | 8/8 | 0.892 |
| middle | 0.75 | 10 partitions, rest 39 | 7/8 | 0.804 |
| middle | 0.70 | 1388 / 599 / rest 13 | 7/8 | 0.844 |
| middle | 0.65 | 1904 / 96 | 8/8 | 0.870 |
| old (15 %) | none | 1 | 7/8 | 0.710 |
| old | 0.80, 0.75 | everything in rest | 7/8 | 0.70 |
| old | 0.70 | rest 1746 + 3 small | 7/8 | 0.685 |
| old | 0.65 | 11 partitions, rest 756 | 7/8 | 0.598 |

Partitioning costs little: a subfamily split over two partitions gives two sets of consensuses
(redundancy, not loss), and two subfamilies sharing a partition are separated by the k-mer tree
inside it. The threshold can therefore be chosen for *size*, not for biology: start at 0.80 and
lower it (or re-cluster the biggest partition at a lower `--id`) until every partition is below
the cap (~20,000 for the k-mer tree, more with `-P`). VSEARCH on 20,000 × 300 bp: 5 s.

## 4. LINEs: the problem is truncation, not length or divergence

Simulated LINE-like family: 6 kb masters, 8 subfamilies that differ by only 0.67 % per tree
branch (40 sites), copies 10 % diverged, 2,000 copies; in the truncated set each copy is a 3′
fragment (exponential, mean 1.5 kb, min 300 bp, 49 full-length).

| ordering | copies | chunk purity |
|---|---|---|
| k-mer Jaccard k = 6 / 8 / 10 | **full-length** | 0.968 / 0.994 / 0.992 |
| k-mer Jaccard k = 5 / 6 / 8 | **5′-truncated** | 0.442 / 0.455 / 0.572 |
| `mafft --retree 0` | 5′-truncated | 0.457 |
| k-mer containment (shared / smaller set) k = 5 / 6 / 8 | 5′-truncated | 0.337 / 0.357 / 0.654 |
| Jaccard on the last 300 / 600 bp only, k = 5–8 | 5′-truncated | 0.55 – 0.62 |
| Jaccard k = 8, 500 longest copies | 5′-truncated | 0.596 |
| Jaccard k = 8, 500 longest copies **trimmed to their common 2.1 kb span** | 5′-truncated | 0.770 |

Full-length 6 kb copies at this divergence order almost perfectly (0.99); the same copies
truncated order at random-plus (0.45). Two effects, both confirmed:

- **Length confounds every whole-sequence distance.** A 400 bp fragment and a 3 kb fragment of
  the same subfamily share few k-mers relative to their union (Jaccard) and the fragment is
  contained in every subfamily (containment). Trimming copies to a common span removes the
  confound (0.60 → 0.77).
- **Information scales with shared span.** Sister subfamilies differ at 1.3 % of sites: 4 sites
  in a 300 bp fragment, 28 in 2.1 kb, 80 in the full element. Short fragments cannot be placed
  by any method; they carry no subfamily signal.

The consensus step has the matching problem: with `-p 0.36` a column is called only where ≥ 36 %
of the chunk has a base, so the consensus of a chunk of fragments ends where coverage thins
(mean consensus length 1.6 kb here). SINEderella's `sine_consensus.sh` already uses the right
rule for this: plurality over the *covering* sequences (gaps excluded from the denominator), plus
a minimum coverage.

### Design consequences for long, truncated elements

1. **Discover from the long copies, assign the fragments.** Take the copies above a length
   quantile (or above a fixed 3′ span), trim them to the common span, run SubFam on those. Then
   assign every fragment to the consensus it matches best (`vsearch --usearch_global --iddef 2`
   or ssearch36, exactly asSINEment's job). Fragments are evidence for a subfamily, not
   discoverers of it.
2. **Compare over the overlap, never over the union.** Ordering distance for mixed-length copies
   must be computed on the shared interval (position-aware k-mers after mapping to a backbone, or
   `iddef 2` alignment identity on a subsample), not whole-sequence Jaccard.
3. **Coverage-relative plurality** (`-p` as a fraction of covering sequences, with `-c MIN`
   minimum coverage) so chunk consensuses extend as far as the data does.
4. **Expect scalar alignment speed above 5 kb** in VSEARCH, and O(L²) per chunk in MAFFT: a
   50 × 6 kb chunk takes seconds, fine; 2,000 such chunks is minutes on a workstation.

## 5. Proposed SubFam-XL

```
copies (10^5 - 10^6, any length)
  │
  ├─ 1. length stratum: copies ≥ S bp (S = a 3′ span the family's long copies share), trimmed to S
  │       short fragments set aside for step 5
  ├─ 2. partition: vsearch --cluster_fast, --id lowered until max partition ≤ cap (§3)
  │       clusters < N pooled into "rest"
  ├─ 3. per partition: SubFam (k-mer tree ≤ 20k, else -P), coverage-relative plurality
  ├─ 4. level 2: SubFam on all consensuses (or per-partition L-INS-i when ≤ ~500 rows)
  └─ 5. assign every copy, fragments included, to its best consensus (usearch_global --iddef 2);
          chunks.tsv (both levels) + assignments = full lineage of every copy
```

Steps 2, 3 and 5 exist today (vsearch, SubFam 1.1, asSINEment); step 1 is a `seqkit`/awk line;
step 4 is SubFam called twice. The only new code is the coverage-relative plurality (`-c`), a
one-line change in the awk consensus. The sparse k-mer neighbour tree (§2, "counter scan") would
later remove the partition step and the 20k cap of the k-mer tree, but it is not needed to get
to a million copies.
