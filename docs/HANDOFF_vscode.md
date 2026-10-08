# Handoff to Claude Code in VS Code

Written 2026-10-08 at the end of a long cloud session. Read this first. Everything marked VERIFIED was read in a source or run; UNVERIFIED is a belief.

## 0. How to work with the owner (Toki / Sergei)
- Short answers. Verify from source text or code, never from memory. Separate VERIFIED from UNVERIFIED.
- Do not flatter. "SINEderella is much better than COSEG" is the owner's hypothesis, not a result.
- The owner corrects mistakes quickly. Check repos before saying something is missing (twice in this session a "missing" thing was in a repo or a different tool).
- When naming something that must be fetched, give a link.
- No model IDs in repo content. Commit trailers per the session system reminder. Do not open a PR unless asked.
- Compute runs on the owner's servers (DRAGEN, KIT: see SINE-discriminator/CLAUDE.md). The cloud container is ephemeral and cannot reach NCBI, UCSC, Dfam, Wayback, publishers; GitHub and apt work.
- Update this file after each exchange.

## 1. Repos and branches (all github.com/Toki-bio)
| repo | branch | state |
|---|---|---|
| SubFam | `claude/practical-dijkstra-a0xv8i` | SubFam 1.2.0 + docs + benchmark/alu_konkel; no PR. Earlier work (v1.1.0) is also on `ccr-85af825c-on1o5u`. |
| sinederella | `subfam-1.2` | SubFam wrapper (see 3). From main d418439. No PR. `step1c-deplete` (a2253adb, deplete loop + SubFam patch) is untouched; main untouched. |
| SINE-discriminator | main (read only) | docs read, code not. |
| coseg (rmhubley) | read only | built and run here. |

## 2. SubFam 1.2.0 (what exists)
`SubFam.sh`: k-mer guide-tree ordering (port of ViewAlign kmer-tree.js), chunks of N, MAFFT per chunk, plurality consensus (fraction -p 0.36), final L-INS-i alignment of consensuses. Options: -n -p -k -t -o -x -r -m -P -a -K -v -h, and new `-c` (coverage-relative plurality for truncated copies).
- `-c` VERIFIED on simulated 5'-truncated copies (consensus lengths 462-1493 bp default vs 513-1499 with -c; full-length input gives 1500 bp for every chunk); default path VERIFIED unchanged against v1.1.0 (consensus sequences identical, headers differ by prefix only). EMBOSS equivalence test (`tests/check_equivalence.sh`) was NOT re-run (EMBOSS not installed here).
- Other docs: README.md, docs/SCALING.md, docs/REALDATA_BENCHMARK.md (benchmark plan, partly superseded by this file), docs/papers/*.md (paper notes, section 5).
- No LICENSE (owner and Nikita to choose; options MIT/BSD-3 permissive, GPL-3 copyleft; ViewAlign is MIT). Name "SubFam": no other software of that name found in web search (GitHub itself could not be searched from the container; check manually).

## 3. SINEderella change (branch subfam-1.2)
`SubFam` is now a wrapper around vendored SubFam 1.2.0 (`tools/vendor/SubFam.sh`) with the old contract: run in the working directory as `SubFam input.fasta [BIN]`; writes `input_NNN.bnk` (chunk sequences), `input_NNN.cons` (consensus, '-' where none), `input.clw` (consensuses concatenated, unaligned), `input.msf`, plus new `input.chunks.tsv`. Old script kept as `tools/vendor/SubFam.old.sh`. Env: THREADS, SUBFAM_ORDER=mafft (old ordering), SUBFAM_PLURALITY, SUBFAM_COVERAGE=1. README and MANUAL 6.1.1 updated.
- VERIFIED: runs in a scratch dir on 316 Alu bodies, default and SUBFAM_ORDER=mafft, files and names as above.
- NOT tested: inside a full SINEderella run (step1, step8a, extract_alignments call `SubFam input.fasta 50`); comparison with the old script's output (needs EMBOSS); behaviour with `-r` strands. Default k-mer ordering changes chunk composition vs old runs (SUBFAM_ORDER=mafft restores it).
- Decision for the owner: merge this branch? And merge `step1c-deplete` (only its plurality change alters results, only for bin sizes other than 50)?

## 4. Findings from this session
### 4.1 The literature (details in docs/papers/)
- Price, Eskin, Pevzner 2004 (full text read): one reference (AluSx), copies missing >5 bases at either end dropped (~480k of >1M), indels excluded from the tests, de novo split from one subfamily to a 60-subfamily scaffold, then 213 (+153 by single mutations). Owner's criticism (one collinear reference, no truncated copies, no indels) CONFIRMED from the Methods. "Only works for tight known families" is not demonstrated by the paper either way. The assignment rule and the age-bias formula are not in the text. Supplement (213 consensuses) received and parsed (213 subfamilies, 476,152 copies, matches Table 2).
- COSEG (README read, run here) is Price's code, maintained by Hubley/Smit/Siegel; needs alignments to one reference; Dfam curators say it is "not appropriate" for truncated-looking sets and recommend visual review; cd-hit pipeline is their alternative for length variants (owner: cd-hit is the wrong tool for diverged SINEs).
- Carey et al. 2020 (preprint, not peer reviewed): >10 % of replicate Alu/L1 copies get different RepeatMasker subfamilies.
- Ahmed 2013, Styles & Brookfield 2009, Kryatova 2017, Konkel 2015, Gardner/MELT 2017: young subfamilies defined by one or two diagnostic changes plus a minimum number of independent copies, no statistical tests; gene conversion and homoplasy documented; nomenclature inconsistent (Konkel: one name with several consensuses, one consensus with several names). Orthology/polymorphism evidence is a separate dimension (owner has SINE_orth_loc/ComPair).
- Correction made: an earlier "10-20 % gene conversion" figure was wrong (it is the share of secondary source genes); fixed in docs/REALDATA_BENCHMARK.md.

### 4.2 Real-data tests (benchmark/alu_konkel/, README there has the tables)
- Data: the 343 GenBank records KT305395-KT305737 (Konkel 2015 Sanger loci). Fetch: see section 7.
- SubFam chunks vs lineage labels (labels = best Price consensus among AluY/Ya5/Yb8/Yb9; similarity-based, NOT independent truth): purity 0.918 (n=20), 0.937 (n=10) vs ~0.53 random.
- COSEG (-k -d, our own alignment conversion, 295 copies) vs SubFam chunks on the same copies: COSEG separates Yb8 from the rest but never Ya5 from AluY (even at -m 5); at equal granularity COSEG scores higher (V 0.62-0.68 vs 0.45-0.50), SubFam -n 20 is very homogeneous (0.93) but over-split by design (completeness 0.36). This does NOT show SINEderella is better than COSEG; it tests only the chunk step, not the peel. Reason COSEG missed Ya5 not investigated.

### 4.3 SINEderella and SINE-discriminator (docs read, not code)
- SINEderella: step 2 = flat-bank 10-cycle `ssearch36` vote (unanimous + 0.45 x 10th-best), alignments not stored; the subfamily call is made by a person on SubFam output (MANUAL 6.1) or the peel in SINE-discriminator (`SUBFAMILY_METHOD.md`: features = (position, character incl. gap), Jaccard blocks, peel; purity 0.937 on Timema vs his curation, ceiling 0.882; `best` split/merge threshold calibrated on two cases). docs/FAMILY_SUBFAMILY_ASSIGNMENT.md states whole-length bitscore is the wrong instrument for subfamilies. No repo cites Price or COSEG. CpG handling is only an open divergence decision (D2); not found in the peel.
- The owner's view: no good instrument exists; the aim is to hand an expert the best recalculated data.

## 4.4 Addendum (later in the same day): the Alu history doc and COSEG source
- The owner's living doc "Alu subfamily definitions and how they are told apart" is a Claude Docs artifact: https://claude.ai/artifact/U5kNLaQtKsZByEy5sr6fRD (read and edited through the Claude Docs connector, not the Artifact tool; doc id db5071ff-dd92-4565-8716-30e2a2e83eec). Added sections 3A (COSEG at source level) and 8A (code-level comparison with the peel and the assignment vote), corrected section 8 (the split/merge rule is in pairsurvey.py, not the peel; 0.80 in code vs 0.85 in docs), updated section 9 (sources not read). Details: docs/papers/coseg_source.md.
- COSEG default p-value: integer division makes the correction constant 1 (verified by running an instrumented copy; no effect on tested data). Not reported to the authors.
- Peel code facts (read, not run on Alu): see coseg_source.md last paragraph.

## 4.5 Decisions of 2026-10-08 (evening)
- Genome hg38. Reference sets for the SubFam paper's real-data claim: Dfam Alu families AND Price's 213 consensuses.
- Paper split: SubFam paper = tool + simulated/scaling results + consensus recovery on hg38 Alu (benchmark/alu_hg38/, protocol and metrics fixed in its README before any run); no adjudication, no subfamily-calling claim, no COSEG claim. SINEderella paper = discovery/assignment, comparison with RepeatMasker (stated Dfam version), COSEG, adjudication of disagreements by independent evidence.
- Words: "adjudicate" (decide which call independent evidence supports), "better supported" (not "right"), "demonstrable claim".
- NOT yet run: nothing in benchmark/alu_hg38 has been run on hg38 (needs the server). recover.py was only script-checked on the Konkel loci.

## 5. Paper notes in docs/papers/
price2004.md (incl. supplement check), carey2020_storer2021.md, ahmed2013.md, styles2009.md (incl. S3-S7 counts: 55/3/11/233/91 records; S1, S2 missing; two files hold fewer records than the paper's tables, unexplained), kryatova2017_konkel2015_gardner2017.md.

## 6. Verdict given on a SubFam manuscript (no draft exists in any repo)
Defensible now: a threshold-free compression tool for large copy sets, simulated benchmarks, 100k scaling run, first coarse real-data result. Not yet supported: calling it "subfamily identification", superiority to COSEG, fine-level recovery on real loci, old (AluS/AluJ) families. Reviewer questions to expect: why Price/COSEG are not cited, what independent truth was used, chunk smoothing of minority diagnostics (documented in the owner's own notes). If a draft exists, ask the owner to upload it.

## 7. Reproducing / data that is NOT in the repo
- Konkel loci (NCBI blocked in the container; the owner's server works with retries):
  ```
  mkdir -p konkel && cd konkel
  for i in $(seq 395 737); do for try in 1 2 3 4 5; do
    curl -sf -m 60 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KT305$i&rettype=gb&retmode=text" -o KT305$i.gb && break; sleep $((try*2)); done; sleep 0.4; done
  ```
  (343 records, ESearch count verified = 343.) The owner already did this once (NCBI.ZIP).
- Price supplement `SUBFAMILIES.FINAL.txt` (publisher file, not committed; owner has it; article page DOI 10.1101/gr.2693004 "Supplemental Material").
- COSEG: `git clone https://github.com/rmhubley/coseg; make`; needs perl. Its `ALU.cons` is the AluSx consensus used for the Konkel run.
- Pipeline: `prep.py GB_DIR PRICE_SUPPLEMENT OUT` (needs biopython) -> `SubFam.sh -n 20 alu_bodies.fa`; `to_coseg.py ALU.cons alu_bodies.fa konkel` -> `runcoseg.pl -k -d -m M -c ALU.cons -s konkel.seqs -i konkel.ins`; `compare.py loci_labels.tsv names295.txt FILE:LABEL ...`.
- Container tools were installed with apt (mafft, gawk) and pip (biopython); mawk is the default awk.

## 8. Not done / open
1. Fine-level Alu truth: the authors' per-locus subfamily table (Konkel supplementary S3) or our own diagnostic-position calls; needed for Yb7a3/Yb8b1/Ya4a1/Yb10/Yb11 (3-16 copies, below SubFam's default chunk of 50).
2. Test the peel (not just chunks) against COSEG and Price's partition on the same copies; full-genome Alu (hg16/hg38) needs the owner's server and is where Price's 480k-copy scale applies. Metrics were fixed after the first look this time; for the next test fix them in advance and do blind manual review of disagreements in MSA-viewer.
3. Why COSEG did not split Ya5 from AluY on the Konkel data (alignment conversion, blacklist, minimum size?).
4. Run the SINEderella branch inside a full run; compare with old output.
5. CpG decision D2; whether the peel needs a CpG rule.
6. LICENSE, name check on GitHub, Zenodo DOI, optional Bioconda recipe; reword the SINEderella README line that said SubFam does "subfamily identification" (done on the branch).
7. Blocked sources still wanted: Wheeler lab data (Carey), alucode tarball (403 on every mirror tried; its code is inside COSEG), Dfam pages, Styles S1/S2 alignments, Konkel S3 table.
8. Coverage-relative plurality is implemented; hierarchical second pass for >100k inputs and a sparse k-mer neighbour tree (>20k) are not.
