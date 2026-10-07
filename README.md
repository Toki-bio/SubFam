# SubFam

**Reduce thousands of related sequences to a small alignment of consensus sequences, without choosing an identity threshold.**

SubFam is meant for datasets made of many copies of the same thing that have diverged into internal groups:

- **transposable elements**: all genomic copies of a family (Alu, L1, SINEs, MITEs, …) → its subfamilies
- **satellite DNA**: thousands of monomers → monomer variants and higher-order-repeat units
- **multigene families**: rDNA, histone, tRNA, olfactory receptor copies → paralogous groups

The output is a short alignment (for example 40 rows for 2,000 copies) where every row is the consensus of a group of closely related copies. The private mutations of individual copies are averaged away, and the variants a group shares stay visible. You can read that alignment by eye, build a tree from it, or use it as a library for annotation.

## How it works

```
 2,000 copies         k-mer tree order        chunks of N        chunk consensus       final alignment
 (unordered)    ─►    (similar copies    ─►   neighbours    ─►   (mafft + plurality ─► (mafft L-INS-i)
                       end up adjacent)       N = 50             consensus)
```

1. **Order** all sequences along a k-mer guide tree. The distance between two copies is the weighted Jaccard distance of their k-mer counts (default k = 6). The tree is built by UPGMA, and at every merge the two subtrees are flipped so that their closest ends meet. No alignment is computed in this step, only an order in which related copies sit next to each other. The method is a port of the "Reorder by similarity" guide tree in [ViewAlign](https://github.com/Toki-bio/MSA-viewer) (`kmer-tree.js`). `-m` uses the MAFFT guide tree instead (`mafft --retree 0 --reorder`, as in SubFam 1.0).
2. **Chunk** that order into consecutive blocks of *N* sequences. Leftover sequences (fewer than *N*) join the last chunk, so nothing is dropped.
3. **Consensus**: align each chunk (MAFFT, in parallel) and call its plurality consensus. A base is called when at least `-p` × chunk size of the sequences agree. The consensus is EMBOSS `cons` reimplemented in awk (EDNAFULL scores, the same tie rules). Its output is byte-identical to `cons`, which is no longer needed.
4. **Align** the chunk consensuses with MAFFT L-INS-i.

### Why chunks and not clusters?

Standard tools (CD-HIT, VSEARCH, MMseqs2) build clusters around an **identity threshold**. For repeat families that threshold rarely works:

- **Young families** (copies about 3% from their source): a 95% threshold splits them into hundreds of tiny clusters, and an 80% threshold lumps all subfamilies together.
- **Old families** (copies about 15% from their source, so two copies are about 25–30% apart): any threshold that separates subfamilies leaves most copies as **singletons**.
- The one threshold that does work is different for every family. You can't know it in advance, and you can't apply a single one to a dataset that mixes young and old families.

SubFam never asks "are these two copies X% identical?". It only asks "which copies are each other's nearest neighbours?" and then averages *N* of them. Averaging works at any age: with 50 copies, a consensus recovers the source sequence even when each copy differs from it at 15% of positions. A chunk boundary may still fall inside a subfamily, or straddle two. You then get a few redundant or mixed rows, but each consensus is still exact for the majority of its chunk.

**Choosing N:** the most abundant subfamily dominates every chunk it fills, so a subfamily needs roughly ≥ *N* copies to get its own row. Lower *N* resolves rarer subfamilies at the cost of more rows and noisier consensuses.

## Installation

Requirements: `bash`, `awk` (POSIX; tested with gawk and mawk), [MAFFT](https://mafft.cbrc.jp/alignment/software/) ≥ 7, and `python3` with `numpy` for the k-mer ordering (not needed with `-m`). EMBOSS is no longer required.

```bash
conda install -c conda-forge -c bioconda mafft numpy
git clone https://github.com/toki-bio/SubFam && cd SubFam
./SubFam.sh -h
```

## Usage

```bash
./SubFam.sh [options] sequences.fasta
```

| option | default | meaning |
|---|---|---|
| `-n INT` | 50 | sequences per chunk |
| `-p FLOAT` | 0.36 | fraction of a chunk that must agree to call a base |
| `-k INT` | 6 | k-mer size of the ordering tree (3–12) |
| `-t INT` | all cores | threads |
| `-o DIR` | `.` | output directory |
| `-x STR` | input name | output prefix |
| `-r` | off | copies may be on both strands: strand-independent (canonical) k-mers for the tree, then each copy is oriented like its neighbour in the order |
| `-m` | off | order with the MAFFT guide tree instead of k-mers |
| `-P` | off | order with MAFFT PartTree (implies `-m`), for inputs too large for an all-against-all matrix |
| `-a` | off | keep no-consensus positions as `N` (by default they are removed) |
| `-K` | off | keep intermediate chunk files |

Outputs:

| file | content |
|---|---|
| `PREFIX.cons.fasta` | one consensus per chunk, unaligned |
| `PREFIX.aln.fasta` | the consensuses aligned (FASTA) |
| `PREFIX.msf` | the same alignment in MSF format |
| `PREFIX.chunks.tsv` | `input id`, `consensus`, `strand`: which copy built which consensus |

`annotate_aligned_by_ssearch_seqkit.sh` labels the resulting consensuses with their best hit in one or more reference libraries (Dfam, RepBase, your own).

## How SubFam compares to the obvious alternatives

Each of these tools can produce "fewer, representative sequences". The commands below show the closest replacement for SubFam that each one offers:

```bash
# VSEARCH: greedy centroid clustering, with its own built-in consensus output
vsearch --cluster_fast copies.fa --id 0.80 --consout reps.fa --uc clusters.uc

# CD-HIT-EST: greedy clustering; representatives are real copies (centroids)
cd-hit-est -i copies.fa -o reps.fa -c 0.90 -n 8
#   ...and to get consensuses you must align every cluster and call cons yourself
#   (benchmark/cluster_consensus.sh does exactly that)

# MMseqs2: cascaded clustering; representatives are real copies
mmseqs easy-cluster copies.fa reps tmp --min-seq-id 0.90 -c 0.8
```

They differ from SubFam in two ways. First, every one of them needs `--id`/`-c` chosen per dataset. Second, CD-HIT and MMseqs2 return a *real copy* as representative, and a real copy carries all of its own private mutations (about 15% for an old family). The benchmark below measures what that costs.

## Benchmark

`benchmark/` contains a simulator and a script that runs every method on the same input:

```bash
benchmark/run_benchmark.sh bench_out 8      # needs mafft, numpy, emboss (reference consensus for the clustering tools), cd-hit, vsearch, mmseqs2
```

**Simulated families.** Each family starts from a 300 bp ancestor and evolves along a tree into 8 subfamilies. Every branch adds 4 diagnostic substitutions, so sister subfamilies differ at 8 positions (2.7%). The subfamilies have 600, 400, 300, 250, 200, 150, 60 and 40 copies (2,000 in total). Each copy then decays independently from its subfamily source (85% substitutions, 15% short indels) by 3% (young), 8% (middle) or 15% (old). Copies are shuffled and renamed, so no method sees the labels.

**Metrics.** *reps*: number of output sequences. *reps≥10*: those built from at least 10 copies. *singletons*: fraction of copies left alone in a cluster. *median id*: identity of the reps≥10 to the nearest true subfamily source. *recovered*: true subfamilies matched at ≥ 99% identity by a representative whose group is mostly that subfamily.

Each cell shows **subfamilies recovered (of 8)** and, in brackets, **how many representatives** the method produced from 2,000 copies. Full metrics are in [`benchmark/results.tsv`](benchmark/results.tsv).

| method | young (3%) | middle (8%) | old (15%) |
|---|---|---|---|
| **SubFam `-n 50`** (no threshold) | **8/8** (40) | **8/8** (40) | **7/8** (40) |
| **SubFam `-n 20`** (no threshold) | **8/8** (100) | **8/8** (100) | **7/8** (100) |
| VSEARCH `--consout --id 0.95` | 8/8 (1,306) | 0/8 (2,000) | 0/8 (2,000) |
| VSEARCH `--consout --id 0.90` | 8/8 (97) | 1/8 (1,961) | 0/8 (2,000) |
| VSEARCH `--consout --id 0.85` | 6/8 (6) | 8/8 (1,133) | 0/8 (2,000) |
| VSEARCH `--consout --id 0.80` | 0/8 (1) | 8/8 (210) | 0/8 (1,988) |
| VSEARCH `--consout --id 0.75` | 0/8 (1) | 7/8 (15) | 3/8 (1,536) |
| VSEARCH `--consout --id 0.70` | 0/8 (1) | 0/8 (3) | 7/8 (571) |
| VSEARCH `--consout --id 0.65` | 0/8 (1) | 0/8 (2) | 2/8 (73) |
| CD-HIT-EST `-c 0.95` + mafft/cons per cluster | 8/8 (650) | 0/8 (2,000) | 0/8 (2,000) |
| CD-HIT-EST `-c 0.90` + mafft/cons per cluster | 4/8 (10) | 7/8 (1,660) | 0/8 (2,000) |
| CD-HIT-EST `-c 0.80` + mafft/cons per cluster | 0/8 (1) | 6/8 (21) | 4/8 (1,559) |
| CD-HIT-EST centroids, `-c` 0.95 / 0.90 / 0.80 | 0/8 | 0/8 | 0/8 |
| MMseqs2 centroids, `--min-seq-id` 0.95 / 0.90 / 0.80 | ≤ 3/8 | 0/8 | 0/8 |

What the table shows:

- **A threshold that works for one family fails for another.** VSEARCH matches SubFam on each family only inside a narrow window: 0.90 for young, 0.80–0.85 for middle, exactly 0.70 for old. One step outside that window, it either splits the copies into singletons or merges all subfamilies into a single row. A real genome mixes families of every age, so no single `--id` can be right.
- **Even inside that window it is less compact.** On the old family, the best VSEARCH setting needs 571 representatives (16% of copies left as singletons) to reach what SubFam gets with 40.
- **Real copies are bad representatives.** CD-HIT and MMseqs2 centroids never recover a subfamily source. They carry their own mutations: median identity to the true source is 94–98% (young) and as low as 86–90% (middle and old). Building a consensus per cluster fixes this, but only at the cluster sizes the threshold happens to produce.
- **SubFam's consensuses are exact.** Median identity to the true source is 100.0% in all three scenarios. The single miss is the smallest subfamily (40 copies, below *N*), which merges into its sister subfamily: see [Limitations](#limitations).
- **The cost is speed.** SubFam takes 7–15 s on 2,000 × 300 bp (4 cores), versus about 1 s for VSEARCH.

## Limitations

- Fixed-size chunks do not follow subfamily boundaries. Expect some redundant rows (an abundant subfamily filling several chunks) and some mixed rows (a chunk straddling two subfamilies). Use `PREFIX.chunks.tsv` to see which copies built each row.
- Subfamilies with fewer than about *N* copies may be absorbed by a neighbour. Lower `-n` if rare variants matter.
- Consensus calling assumes the copies are roughly collinear. Strongly truncated or rearranged copies (for example 5′-truncated LINEs) should be trimmed or filtered first.
- The ordering step builds an all-against-all distance matrix. Use `-P` (PartTree) for very large inputs.

## Authors

Original idea: Toki. Parallel chunk alignment: Nikita.

## License

Not yet chosen. Until a LICENSE file is added, all rights are reserved by the authors.
