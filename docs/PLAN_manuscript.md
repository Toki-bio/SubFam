# Plan: the SubFam manuscript (revised 2026-10-10)

Inputs: the discussion of 2026-10-10, `docs/research/B2_literature.md`,
`docs/research/RELATED_WORK.md`, `docs/research/CODE_REVIEW.md`. Nothing below has been started.

## 0. Decisions already taken

- Scope: two objectives. (1) A visually readable presentation of a variety of similar sequences,
  as similarity subgroups. (2) An alignment proxy that lets downstream tools handle thousands of
  repetitive sequences. One sentence of limits beside them: no identity threshold, a chunk is not a
  subfamily, orphans and truncated copies are not resolved.
- Framing sentence (author, 2026-10-10): SubFam does not claim start-to-end work; it bridges the
  gap between the gold standard of manual curation and a researcher's computational analysis by
  organising the pile of chaotic raw sequences into an ordered, readable form.
- Journal: BMC Bioinformatics (Software) or Mobile DNA. Nature Communications only if the B2 example
  yields a biological finding worth its own paper; as a tool paper with three examples it would be
  asked for the broad benchmarking this plan avoids.
- Reproducibility: every example is rebuilt by `examples/rebuild.sh` from named inputs and seeds;
  every benchmark directory has a README with the exact commands.
- Three working examples: B2 (main figure), L1 (long, truncated copies), Alu (crowded family,
  supplement). A satellite example is deferred until B2 is on screen.
- Simulations stay, demoted to one paragraph + one supplementary figure: they test one property
  (CpG decay / old families) that real data cannot.
- SINEderella is not exposed. The peel is removed from the paper's L1 example.
- Downstream integration: a short theoretical section; no new method, no new claim. One cheap concrete
  case (a tree of B2 chunk consensuses) kept in reserve for reviewers.
- Benchmarking kept minimal: one metric, three baselines (section 5).

## 1. Manuscript skeleton (BMC Bioinformatics, Software)

1. Background — the multitude problem; why thresholds fail across ages; what a proxy alignment is.
2. Implementation — order, chunk, consensus, align; parameters `-n -p -k -c`; what 0.36 is (18/50) and
   what it does; complexity and limits (dense n² matrix: fine to ~2×10⁴, `-P` beyond).
3. Results —
   3.1 B2 (figure 1: the grouped alignment, heads + chunk rows);
   3.2 L1PA2–8 (figure 2: ragged 5′ ends, classes at ORF2 resolution, the 3′-end control);
   3.3 Alu (supplementary figure; Price consensuses beside their groups);
   3.4 Baselines (one table: vsearch sweep, single consensus, [COSEG]);
   3.5 Simulation (supplementary: CpG decay, hard-old).
4. Discussion — what the proxy alignment can feed (section 7); what SubFam does not do.
5. Availability — GitHub, MIT, Zenodo DOI; `examples/` with viewer links.

## 2. Example 1: mouse B2 (main figure)

### What the literature gives us (details in `docs/research/B2_literature.md`)

Two classifications of the same family exist, and they disagree in a way that is useful to us:

| Scheme | Units | Basis | Our reading |
|---|---|---|---|
| Repbase/Dfam/rmsk | B2_Mm1a, B2_Mm1t, B2_Mm2, B3, B3A | single-nucleotide diagnostics; Dfam cites Bains 1989 (71 seqs) + MGSC 2002, no subfamily paper | weak provenance; Mm1a/Mm1t differ by 2 nt |
| Vassetzky…Kosushkin, Kramerov 2021 | B2a, B2b, B2c, B2d, B2e | long indels (11-bp central deletion; 28-bp and 9-bp 3′ insertions; 27-bp region before the tail), 92,352 full-length mm10 copies | strong, visible by eye; B2a = Mm1a+Mm1t+Mm2(+Mm1o); B2d = B3; the paper calls Mm1/Mm2 "extreme cases in a continuum" |

So B2 is not merely "few subfamilies": it is a family where an indel-based classification (ours)
and a point-mutation classification (Repbase) coexist. SubFam should reproduce the indel groups
trivially — chunk consensuses keep shared indels — and should show the Mm1/Mm2 continuum as a
gradient across chunk rows rather than as two blocks. That is the figure: both schemes as head
rows, SubFam rows between them, the reader sees which one the data support. SubFam claims nothing;
the 2021 paper already made the claim.

### Data

- Copies: UCSC mm39 rmsk, repFamily = B2 (all five repNames), one chromosome (chr19: whole-genome
  total is 369,746, chr19 ≈ 2–3 % → ~8–9k), keep copies ≥ 170 bp (the 2021 paper used ≥ 90 % length).
  Subsample ~3,000 so that `-n 50` gives ~60 rows — one screen. Files already checked: `rmsk.txt.gz`
  (142 MB) and `chr19.fa.gz` (19 MB) answer HTTP 200. `benchmark/line_real/extract_l1.py` already does
  this cut given a regex.
- Head rows, scheme 1: Dfam consensuses DF000001735–1739 via the API (verified working).
- Head rows, scheme 2: **the B2a–e consensuses of Vassetzky et al. 2021 (Fig. 1A). Needed from the
  author.** Also, if available, the per-copy classification (which of the 92,352 copies is a/b/c/d/e)
  or the "domestic script", so that copies can be labelled by the indel scheme independently of us.
  Fallback if not supplied: label copies by the diagnostic indels ourselves (11-bp deletion, 28-bp /
  9-bp insertions) after aligning to the B2 consensus — workable but it re-derives what the paper did.
- Truth layout as in `examples/`: `TRUE_<group>_consensus_of_<group>_copies` from the real member
  copies, `PUBLISHED_<name>` for Dfam and for the 2021 consensuses.

### Expected result and what would count as failure

- Indel groups (a+b vs c vs d vs e) form clean blocks of chunk rows: success criterion identity ≥ 0.98
  of each chunk row to its own-group head, as in the Alu example.
- Within B2a, rows labelled Mm1a/Mm1t/Mm2 by rmsk majority do not form clean blocks; the diagnostic
  positions vary row by row. That is the expected outcome and supports the 2021 reading. If instead
  Mm1 and Mm2 form clean blocks, the figure still works, and the text says so.
- Runtime: ~3,000 × 190 bp, seconds.

## 3. Example 2: L1PA2–8 (long, truncated copies)

- Replace the peel-based grouping by rmsk labels only: head row = consensus of the real PA-k copies
  (≥ 60 sampled, as now), SubFam rows labelled by majority rmsk class. `examples/build_examples.py`
  already does this given a truth file; `build_groups.py` (peel-based) is no longer used for the paper.
- The resolution observation stays: chunks mix PA2+PA3, PA4+PA5(+6), PA7+PA8, matching Smit 1995
  Table 2 ORF2 classes (`benchmark/line_real/README.md` §2).
- Control: run SubFam on the last 900 bp of the same 2,222 copies (the region the PA names are defined
  on). If chunks then follow the finer PA classes, the resolution explanation is demonstrated.
- Check for length mixing: purity of chunks binned by copy length. If short copies group by length
  rather than class, say so and show that `-c` still gives correct consensuses of the covered region.
- Truncated copies: keep them (that is L1), with `-c`. The ragged 5′ edge is part of the figure.

## 4. Example 3: Alu (supplement)

- Keep `examples/konkel_alu_n20.aln.fasta` as is (Price consensuses as `PUBLISHED_`, median identity
  1.000). Text: the family is crowded, SubFam makes no subfamily claim, the published consensuses sit
  beside their groups. The hg38 Alu run on the owner's servers remains optional.

## 5. Baselines (one table)

Metric everywhere: identity of each SubFam row to the consensus of its own group's real copies (what a
method can reach), plus number of rows.

1. **vsearch threshold sweep** (`--cluster_fast --id 0.80/0.90/0.95/0.98 --consout`) on B2, L1, Alu:
   clusters and singletons per threshold, and how many group heads are recovered by a centroid at
   ≥ 0.98. Shows that no single threshold serves a young and an old family. Easy.
2. **Single family-wide consensus** (EMBOSS-style plurality of all copies): identity of each copy to it
   vs to the nearest SubFam row. Shows what one consensus hides. Easy.
3. **COSEG** on Alu only (its own domain; needs cross_match alignments to a seed). Medium; decide after
   1–2 are done. If skipped, a paragraph: COSEG asserts subfamilies by co-segregation tests and needs
   full-length copies; SubFam asserts nothing and takes fragments.
4. Paragraphs, no runs: PASTA/MAGUS (divide-and-conquer aligners merge subsets into n rows; SubFam
   collapses each subset to one row); TEtrimmer 2025 (tree + DBSCAN split of one family's MSA — the
   nearest published idea, needs a genome and distance parameters); Swarm/ASV tools (error models, not
   mutation); representative-subset pickers (Treemmer/PARNAS pick copies, SubFam averages them).

## 6. Code before any new run (`docs/research/CODE_REVIEW.md`)

Must-fix, then re-run all examples, bump to 1.2.1:
1. mawk crash in the chunk splitter (`int(i++ / n)` → `int(i / n); i++`).
2. Floating-point ceil in the plurality threshold (two places in `SubFam.sh`, one in
   `examples/build_examples.py`). Default `-n 50 -p 0.36` is exact; `-n 20` and all `-c` columns are
   affected, so the L1 example must be re-run after the fix.
3. Prefix from basename before stripping the extension.
4. Glob-safe chunk file listing; die if no consensus was produced.
5. numpy ≥ 1.20 check.
Should-fix in the same release: balanced chunk sizes (last chunk currently N…2N−1), `--threadit 0`,
drop the no-op `--ep 0.123`, chunk name in failure messages, gzip detection.
README: document lowercase = ≤ 50 % support, gap veto, the origin of 0.36, actual chunk-size range,
numpy version, remove the mawk claim or fix it.

## 7. Downstream integration (discussion section, theoretical)

Each item one sentence of what the proxy alignment feeds and why; no runs:
- profile HMM (`hmmbuild` → `nhmmer`) for annotation of a family without a curated library;
- subfamily consensus curation (Dfam/RepeatModeler Refiner start from a short alignment);
- group-wise divergence landscapes without a prior library;
- diagnostic-position discovery (columns bimodal across rows — Smit's criterion, by machine);
- recombinant/mosaic copies (a row whose affinity switches along its length);
- satellite monomer variant maps; probe design from conserved vs variable columns;
- feeding tools that choke on thousands of rows (IQ-TREE, Jalview, Gblocks, ViewAlign).
Reserve: a tree of the B2 chunk consensuses (one IQ-TREE run) if a reviewer asks for a concrete case.
Past publications to cite here: to be named by the author.

## 8. Order of work and what is needed from the author

Order: 6 (code fixes, re-run existing examples) → 2 (B2) → 3 (L1 regroup + 3′ control) → 5.1–5.2
(baselines) → manuscript text → 5.3 (COSEG) only if time.

Needed from the author:
1. B2a–e consensus sequences (Vassetzky et al. 2021, Fig. 1A) and, if available, the per-copy a–e
   classification or the script that makes it.
2. Decision: include truncated L1 copies (recommended) or full-length only.
3. Decision: run COSEG or write the paragraph.
4. The past publications to cite in section 7.
