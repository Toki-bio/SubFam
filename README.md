# SubFam

**Reduce thousands of related sequences to a small alignment of consensus sequences, without choosing an identity threshold.**

SubFam is meant for datasets made of many copies of the same thing that have diverged into internal groups:

- **transposable elements**: all genomic copies of a family (Alu, L1, SINEs, MITEs, …) → its subfamilies
- **satellite DNA**: thousands of monomers → monomer variants and higher-order-repeat units
- **multigene families**: rDNA, histone, tRNA, olfactory receptor copies → paralogous groups

The output is a short alignment (for example 40 rows for 2,000 copies) where every row is the consensus of a group of closely related copies. The private mutations of individual copies are averaged away, and the variants a group shares stay visible. You can read that alignment by eye, build a tree from it, or use it as a library for annotation.

## How it works

```
 2,000 copies         guide-tree order        chunks of N        chunk consensus       final alignment
 (unordered)    ─►    (similar copies    ─►   neighbours    ─►   (mafft + EMBOSS  ─►   (mafft L-INS-i)
                       end up adjacent)       N = 50             cons, plurality)
```

1. **Order** all sequences along the MAFFT guide tree (`mafft --retree 0 --reorder`). This step computes no alignment, only an order in which related copies sit next to each other.
2. **Chunk** that order into consecutive blocks of *N* sequences. Leftover sequences (fewer than *N*) join the last chunk, so nothing is dropped.
3. **Consensus**: align each chunk (in parallel) and call its plurality consensus. A base is called when at least `-p` × chunk size of the sequences agree.
4. **Align** the chunk consensuses with MAFFT L-INS-i.

### Why chunks and not clusters?

Standard tools (CD-HIT, VSEARCH, MMseqs2) build clusters around an **identity threshold**. For repeat families that threshold rarely works:

- **Young families** (copies about 3% from their source): a 95% threshold splits them into hundreds of tiny clusters, and an 80% threshold lumps all subfamilies together.
- **Old families** (copies about 15% from their source, so two copies are about 25–30% apart): any threshold that separates subfamilies leaves most copies as **singletons**.
- The one threshold that does work is different for every family. You can't know it in advance, and you can't apply a single one to a dataset that mixes young and old families.

SubFam never asks "are these two copies X% identical?". It only asks "which copies are each other's nearest neighbours?" and then averages *N* of them. Averaging works at any age: with 50 copies, a consensus recovers the source sequence even when each copy differs from it at 15% of positions. A chunk boundary may still fall inside a subfamily, or straddle two. You then get a few redundant or mixed rows, but each consensus is still exact for the majority of its chunk.

**Choosing N:** the most abundant subfamily dominates every chunk it fills, so a subfamily needs roughly ≥ *N* copies to get its own row. Lower *N* resolves rarer subfamilies at the cost of more rows and noisier consensuses.

## Installation

Requirements: `bash`, `awk`, [MAFFT](https://mafft.cbrc.jp/alignment/software/) ≥ 7, and [EMBOSS](https://emboss.sourceforge.net/) (`cons`, `seqret`).

```bash
conda install -c conda-forge -c bioconda mafft emboss
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
| `-t INT` | all cores | threads |
| `-o DIR` | `.` | output directory |
| `-x STR` | input name | output prefix |
| `-r` | off | copies may be on both strands (MAFFT `--adjustdirection`) |
| `-P` | off | MAFFT PartTree ordering, for more than about 10,000 sequences |
| `-a` | off | keep no-consensus positions as `N` (by default they are removed) |
| `-k` | off | keep intermediate chunk files |

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
benchmark/run_benchmark.sh bench_out 8      # needs mafft, emboss, cd-hit, vsearch, mmseqs2
```

**Simulated families.** Each family starts from a 300 bp ancestor and evolves along a tree into 8 subfamilies. Every branch adds 4 diagnostic substitutions, so sister subfamilies differ at 8 positions (2.7%). The subfamilies have 600, 400, 300, 250, 200, 150, 60 and 40 copies (2,000 in total). Each copy then decays independently from its subfamily source (85% substitutions, 15% short indels) by 3% (young), 8% (middle) or 15% (old). Copies are shuffled and renamed, so no method sees the labels.

**Metrics.** *reps*: number of output sequences. *reps≥10*: those built from at least 10 copies. *singletons*: fraction of copies left alone in a cluster. *median id*: identity of the reps≥10 to the nearest true subfamily source. *recovered*: true subfamilies matched at ≥ 99% identity by a representative whose group is mostly that subfamily.

_Results table: being generated, see `benchmark/run_benchmark.sh`._

## Limitations

- Fixed-size chunks do not follow subfamily boundaries. Expect some redundant rows (an abundant subfamily filling several chunks) and some mixed rows (a chunk straddling two subfamilies). Use `PREFIX.chunks.tsv` to see which copies built each row.
- Subfamilies with fewer than about *N* copies may be absorbed by a neighbour. Lower `-n` if rare variants matter.
- Consensus calling assumes the copies are roughly collinear. Strongly truncated or rearranged copies (for example 5′-truncated LINEs) should be trimmed or filtered first.
- The ordering step builds an all-against-all distance matrix. Use `-P` (PartTree) for very large inputs.

## Authors

Original idea: Toki. Parallel chunk alignment: Nikita.

## License

Not yet chosen. Until a LICENSE file is added, all rights are reserved by the authors.
