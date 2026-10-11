# hg38 Alu: consensus recovery by SubFam (for the SubFam paper)

Question: do SubFam's chunk consensuses, made from random hg38 Alu copies with no identity threshold and no reference,
contain the known Alu consensuses? Reference sets: **Dfam Alu families** and **Price et al. 2004's 213 consensuses**.
This is the only real-data claim the SubFam paper needs. It does not say which annotation is correct, and it does not
call subfamilies (that is the SINEderella paper, below).

## Fixed before any run (do not change after looking at results)
- **Copies:** full-length Alu from the UCSC hg38 `rmsk` table, repFamily "Alu", length 250-330 bp (`extract_alu.sh`).
  Date and the `rmsk` library version are recorded with the run.
- **Samples:** 30,000 random copies, `seqkit sample -s SEED`, SEED = 1, 2, 3 (three independent samples; report the spread).
- **SubFam:** `-n 20` and `-n 50`, k-mer ordering, other options default (SubFam 1.2.0).
- **Recovery metric (`recover.py`):** local alignment of the reference against every chunk consensus (both strands);
  coverage >= 0.90 of the reference required; identity = identical columns / alignment columns. A reference is recovered at
  T = 0.95, 0.98, 0.99 if its best chunk consensus reaches T. Report all three, plus the median best identity.
- **Eligibility strata:** a reference with few copies cannot be recovered from a 30,000 sample. Expected copies in the sample =
  copies of that label in the full set x 30,000 / total. Report recovery separately for expected >= 50, 10-50, < 10.
  Counts: Dfam families from the `rmsk` repName counts; Price's from the supplement (copies per subfamily; his data are
  hg16 Build 34, ours hg38, so his counts are only a guide).
- **Reverse view:** for each chunk consensus, its nearest reference and identity (redundant rows: several chunks near one
  reference; novel rows: none above 0.95).
- **Baselines, same metric, same samples:** (a) VSEARCH centroids at --id 0.90, 0.95, 0.98 (`vsearch --cluster_fast`),
  (b) SubFam `-m` (MAFFT guide-tree ordering), (c) SubFam with chunks of the same size but random order if a random-order
  option is added (not available in 1.2.0, so left out unless added). Compare recovery against the NUMBER of output
  sequences, not at one setting.
- **Not claimed:** that any reference is "right", that chunks are subfamilies, or anything about COSEG.

## Steps (on the server; this container cannot reach UCSC or Dfam)
1. `extract_alu.sh hg38.fa rmsk.txt.gz hg38_alu_full.fa`
2. Dfam: export all curated `Alu*` consensuses for Homo sapiens as FASTA (Dfam download page or its `famdb.py`; check the
   current command and record the Dfam release) -> `dfam_alu.fa`.
3. `python3 prep_price.py SUBFAMILIES.FINAL.txt > price213.fa` (the authors' supplement; not stored in this repo).
4. `./run_subfam.sh hg38_alu_full.fa 30000 1 run_s1 dfam_alu.fa price213.fa` (then seeds 2 and 3).
5. Baselines (a), (b) and score them with `recover.py` the same way.

## Script check (not a result)
`recover.py` was run on the 31 chunk consensuses of the 316 Konkel Y-lineage Alu bodies against Price's 213: 105 / 213 references
reach 0.95, 52 reach 0.98, 37 reach 0.99, and AluY (0.997), AluYa5 (1.000) and AluYb8 (1.000) are among the recovered. This only shows the
script works: those copies are Y-lineage only (no AluJ/AluS can be recovered) and were labelled by similarity to the same consensuses.

## The second paper (SINEderella): adjudication
Compare RepeatMasker (stated Dfam version) with SINEderella per copy on hg38. For each disagreement collect independent
evidence in a fixed order: diagnostic positions read from the alignment; presence/absence polymorphism (1000 Genomes);
chimpanzee orthology; divergence and CpG consistency; a blind manual review (method labels hidden). Report agreement rate and,
for disagreements, how often each side is supported, including where RepeatMasker wins. Decide metrics and order first.

## Run of 2026-10-11: hg38 chr1, 30,000 copies (`run_chr1_30k.sh`, `extract_alu.py`, `score_runs.py`; outputs in `results/`)
Deviation from the protocol above: one chromosome (chr1: 76,301 full-length Alu copies of 250–330 bp by `extract_alu.py`, a Python counterpart of
`extract_alu.sh`) and one sample (30,000, seed 1, `../sine_b2/subsample.py`), because the whole-genome run and the seeds 2–3 are for the author's servers.
References: the 42 curated Dfam 4.0 human Alu families (`results/dfam_alu.fa`, Dfam API, release 4.0 of May 2026) and the 213 Price et al. 2004
consensuses (from the supplement, kept outside the repository). SubFam 1.4.1 (`-n 50`, `-n 20`, `-n 50 -P`), 8 threads on 4 cores; the `-n 20` run used
the pre-1.4.2 final alignment (34 min, of which most in the L-INS-i refinement of 1,500 rows).

| run | wall | peak RSS | rows (distinct) | purity vs rmsk name (45 names) | Dfam, expected ≥ 50 copies (n = 26): at 0.95 / 0.98 / 0.99 | Price, expected ≥ 50 (n = 74) | Price, 10–50 (n = 90) | Price, < 10 (n = 49) |
|---|---|---|---|---|---|---|---|---|
| `-n 50` | 12 min | 4.3 GB | 600 (584) | 0.472 | 25 / 21 / 17 | 74 / 61 / 37 | 88 / 61 / 37 | 45 / 12 / 5 |
| `-n 50 -P` | 6.7 min | 0.8 GB | 600 (588) | 0.400 | 22 / 20 / 12 | 73 / 57 / 29 | 86 / 57 / 29 | 42 / 10 / 3 |
| `-n 20` | 34 min | 4.3 GB | 1,500 | – | – | – | – | – |

Also: Dfam families with 10–50 expected copies 6 / 5 / 2 (of 6) and with < 10, 9 / 5 / 2 (of 10) with the k-mer order; 463 of the 600 rows are within 0.95 of a Dfam family and 506 within 0.95 of a Price consensus (k-mer order).
Reading: the rows of a 30,000-copy sample contain every Price subfamily that is expected to have 50 or more copies in the sample (74 of 74 at ≥ 0.95,
61 at ≥ 0.98) and 25 of the 26 such Dfam families; recovery falls with the expected copy number, as a chunk needs about N copies of a group.
Purity against the 45 RepeatMasker names is low because many names differ at a few positions (AluSx/AluSx1/AluSx3, AluY/AluYc...), which is the
reason the manuscript uses B2 as the stress case and Alu as one sentence. The `-n 20` run was not scored (its consensuses are in `results/`).
