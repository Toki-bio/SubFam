# Real L1 copies (hg38 chr4, L1PA2-L1PA8): SubFam, peel into groups, group consensuses

Run 2026-10-08 in the cloud session (4 cores). Data: UCSC hg38 `rmsk.txt.gz` and `chr4.fa.gz` (hgdownload.soe.ucsc.edu). Not a benchmark with a truth; a worked example on real copies.

1. `extract_l1.py rmsk.txt.gz chr4.fa.gz l1.fa 1000`: every RepeatMasker L1PA2-L1PA8 copy on chr4 of at least 1 kb, without N, on the element strand. 2,222 copies (L1PA2 179, PA3 428, PA4 393, PA5 372, PA6 220, PA7 483, PA8 147); mean length about 3.3 kb.
2. `SubFam.sh -n 20 -c -t 4 -x l1 l1.fa`: 111 chunk consensuses, 5 min 2 s (`subfam.log`).
3. `peel_features.py` of the SINE-discriminator repository on the chunk consensuses (default settings, `peel.log`): 7 groups over 3 rounds, 105 of 111 chunks, 6 left over (G0).
4. `build_groups.py`: for every group, the consensus of its ORIGINAL copies (a random sample of at most 60, coverage-relative plurality 36 %, MAFFT `--retree 2`), not a consensus of chunk consensuses; then one MAFFT alignment of group heads and SubFam rows, grouped. Output: `examples/real_l1pa_chr4_n20c.aln.fasta`.

Row names: `CONS_G<k>_<sampled>of<all>copies` = group head; `<chunk>|G<k>|<rmsk majority>|p<purity>` = SubFam chunk consensus, the RepeatMasker subfamily that most of its copies carry, and the fraction of its copies carrying it.

## What it shows (read before quoting)
- 105 SubFam rows against the head of their own group: median identity 0.991 (min 0.812, max 1.000) over the shared span, median span 2.6 kb. For 30 of the 105 rows another group's head is marginally closer over the row's span: the heads are not far apart.
- The groups do not reproduce the RepeatMasker labels: chunk purity against rmsk is 0.69 on average; groups G1, G4, G5 mix 3-4 RepeatMasker subfamilies (e.g. G4: PA4 286, PA5 267, PA6 95, PA3 33 of 700 copies); G3 is mostly PA7/PA8, G2 and G7 mostly PA2/PA3. L1PA2-8 are age classes along one lineage, and the rmsk label is itself the best score to a consensus; there is no independent truth here.
- Peel group sizes are uneven (35 chunks vs 5) and G4 was peeled on only 4 features, so it may be several lineages.
- The sample of at most 60 copies per head means G4's head rests on 60 of 700 copies.
Neither "SubFam recovers the L1PA subfamilies" nor the opposite follows from this. It shows what SubFam's output and the peel give on real L1 copies, for the reader to inspect.
