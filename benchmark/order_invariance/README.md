# Why the k-mer order depends on the input order (2026-11-02, DRAGEN, 1-2 threads, SubFam fca3b01 ordering code)

3,000 chr19 B2 copies (make_b2.sh sample, seed 1) and ten shuffles of them (seeds 101-110).

`diag_ties.sh` (instrumented copy of kmer_order.c; its order is identical to the release code): among the 2,999 UPGMA merges the global minimum distance was tied in 144 (and one merge had a tie inside its own row); the ties are broken by position in the file. Against the clean run the ten shuffled trees share 79-96 % of their clades; 81-96 % of neighbouring pairs of the leaf order are shared.

`tb_patch.py` + `diag_tiebreak.sh`: the same code with ties broken by the sequence content (FNV hash of the sequence ranks clusters; tied end pairs of the closest-ends flip are ranked the same way). All ten shuffled trees are identical to the clean tree (clades shared 1.0000); the leaf orders differ at one neighbour pair (0.999; two copies in the sample have the same sequence).

`e2e_tiebreak.sh`: SubFam with that ordering on the clean input and three shuffles: copy pairs kept in the same chunk 1.000, 60/60 original rows recurring at identity 1.000 in all three.

So the order dependence is the tie-breaking, not the closest-ends rule. The release (1.5.1) keeps the file-position tie-break, which is what ViewAlign's kmer-tree.js does (tests/check_equivalence.sh); the content-based tie-break is not in it, because it changes every order and so every number of the manuscript.
