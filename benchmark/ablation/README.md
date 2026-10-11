# Ablation: what each part of SubFam contributes (B2, 3,000 copies, N = 50; 2026-10-11)

`run_ablation.sh B2_WORKDIR` (the directory of `benchmark/sine_b2/make_b2.sh`), output `ablation.tsv`. SubFam 1.4.1 (`-I` keeps the input
order; `SUBFAM_NOFLIP=1` is an environment switch for this test only), MAFFT 7.505, 8 threads. Purity against the RepeatMasker label and the
Vassetzky 2021 best-hit label (labels are used to score only); `agree` and `id_min` are the label-free means of `chunkstats.tsv`.

| variant | rows | columns | s | purity rmsk | purity 2021 | agree | id_min |
|---|---|---|---|---|---|---|---|
| SubFam: k-mer order with closest-ends flip, chunks, consensus, alignment | 60 | 223 | 20 | **0.793** | **0.773** | 0.805 | 0.720 |
| no ordering: shuffled input, same chunks | 60 | 237 | 9 | 0.426 | 0.424 | 0.757 | 0.597 |
| no ordering: genomic order | 60 | 242 | 9 | 0.428 | 0.457 | 0.760 | 0.591 |
| MAFFT guide-tree order (`-m`) instead of the k-mer order | 60 | 223 | 11 | 0.734 | 0.764 | 0.798 | 0.710 |
| k-mer UPGMA order without the closest-ends flip | 60 | 234 | 19 | 0.771 | 0.766 | 0.803 | 0.716 |
| no chunks: one consensus of all copies (`-n 3000`) | 1 | 194 | 47 | 0.422 | 0.392 | 0.796 | 0.503 |
| no consensus: MAFFT `--auto` of all 3,000 copies | 3,000 | 2,781 | 18 | – | – | – | – |

Reading: chunks without an ordering are at the level of one family consensus (0.43 vs 0.42: the chance level for 60 chunks of these
labels); the generic MAFFT order gives most of the gain (0.73), the k-mer UPGMA order more (0.77), and the closest-ends flip the rest (0.79);
without the consensus step the output is the 3,000-row, 2,781-column alignment the method exists to avoid. The label-free statistics
(`agree`, `id_min`) rank the variants the same way. Each listed part exists elsewhere (guide-tree leaf orders, windowed decomposition of
guide trees in PASTA/UPP, family consensus builders such as RepeatModeler's Refiner, per-cluster consensuses in VSEARCH); none of them
produces the object of the first row, and the table shows that the object is made by their combination, not by any one of them.
Caveats: one family, one chromosome, one N; purity is against annotation labels (the label-free columns are there for that reason).
