# Handoff for Claude Code in VS Code (living file; newest facts first in each section)

Last updated 2026-10-08. Keep it short; update after each exchange.

## How to work with the owner
- Short answers. Verify from the source text/code, never from memory or general words.
- Separate VERIFIED (source read) from UNVERIFIED. Do not flatter; do not claim "unique/better" without a test.
- The owner corrects mistakes quickly: check repos before asserting that something is missing.
- Do not put model IDs in repo content. Commit trailers per the session system reminder.
- Compute runs on the servers (DRAGEN, KIT: see SINE-discriminator/CLAUDE.md), not on a laptop.

## Repos (all Toki-bio, GitHub)
- SubFam: branch `claude/practical-dijkstra-a0xv8i` (v1.2.0 = v1.1.0 + `-c`, docs/papers/price2004.md, docs/REALDATA_BENCHMARK.md, this file). Earlier work: `ccr-85af825c-on1o5u`. No PR.
- SINEderella (repo name `sinederella`): branch `step1c-deplete` (a2253adb), main untouched (d418439f).
- SINE-discriminator: read (docs only, not code). Subfamily method: SUBFAMILY_METHOD.md (peel on shared (position, character incl. gap) features).
- COSEG (rmhubley/coseg), read-only reference.
- The cloud container is ephemeral: re-clone, do not expect old paths.

## What was verified this session
- Price et al. 2004 (full text read, notes in docs/papers/price2004.md): one reference (AluSx), copies missing >5 bases at either end dropped (~480k of >1M), indels excluded from the tests, scaffold 60 -> 213 subfamilies, assignment rule and age-bias formula NOT given in the text.
- COSEG README (read): "derived from … programs … written by Alkes Price"; extended to 3 co-segregating mutations; alternative p-value model 2008 (Siegel); "input sequences must be full length alignments to a *single* reference"; `-k` approximates Price's original p-value; repo contains Price's sample data ALU.seqs/ALU.ins/ALU.cons (4,000 copies).
- COSEG built and run here on that sample (`runcoseg.pl -k -d -m 50 -c ALU.cons -s ALU.seqs -i ALU.ins`): scaffold 3 subfamilies (2096/1000/904), 4 overall. NOT yet compared with Price's tables.
- SINEderella: step 2 = flat bank, 10 `ssearch36 -m 8` cycles, unanimous vote + 0.45 threshold; alignments are not stored (m8 has no per-column states). Manual §6.1 = SubFam + human review. docs/FAMILY_SUBFAMILY_ASSIGNMENT.md states whole-length bitscore is the wrong instrument for subfamilies.
- CpG: handled only as a divergence question (decision D2 pending; tools/cpg_div/, docs/CPG_DIVERGENCE_TEST.md, SINE-discriminator/cpg_divergence_review/). Not found in the pipeline or in the peel code/docs. (Owner remembers it as "already in the working": check where.)
- No repo (SubFam, SINEderella, SINE-discriminator) cites Price or COSEG.

## Owner's positions (NOT verified claims)
- No good instrument for subfamilies exists; SINEderella's aim is to give the expert the best recalculated data, the expert looks anyway.
- The peel is meant to be smart and unique but is "not working well enough" yet.
- Believes SINEderella is much better than COSEG. TREAT AS A HYPOTHESIS TO TEST.

## Plan: careful comparison with COSEG (draft, to agree with the owner)
1. Reproduce: COSEG `-k` on the Price sample vs Price Table 1/2 (Ya5 at 0.7%, linkage 12x).
2. Pools: (a) Price/COSEG sample (4,000 Alu); (b) larger Alu from hg16/hg38 (needs server); (c) a curated SINE set, e.g. Timema (13 curated subfamilies; COSEG needs copies projected to one reference, truncated ones dropped: shows the cost).
3. Methods: COSEG (-k and default model, 2 and 3 sites) vs SubFam+peel (+ SINEderella assignment). No tuning on the test pools.
4. Metrics fixed BEFORE running: agreement with curated labels (ARI, per-subfamily recovery incl. small nested ones), fraction of copies classified, runtime.
5. Blind manual review of disagreements in MSA-viewer. Report where COSEG wins as well.

## Blocked from the cloud container (owner offered to fetch; upload PDFs/text)
- Dfam curation guidelines, PMC9191830 (COSEG guidance and the cd-hit-based subfamily pipeline; to be read and cited, even if cd-hit is the wrong tool).
- Research Square rs-86308 (Alu subfamily label reproducibility: >10% of replicates labelled differently, per a search snippet only).
- http://www.cse.ucsd.edu/~ppevzner/download/alucode.tar.gz (403), Price supplement (213 consensuses), Dfam AluSx consensus, UCSC hg16 data.
- GitHub works; apt works (mafft, gawk installed).

## Open items
- LICENSE (owner + Nikita), GitHub name-clash check for "SubFam", merge decision for step1c-deplete (only the plurality change alters results), CpG decision D2.
- SubFam README/SINEderella README wording: SubFam compresses copies; subfamily calls are made by a person.
- Remaining papers to read one by one: Ahmed 2013 (Yb), Styles & Brookfield 2009, MELT, 1000G Alu subfamilies, AluS polymorphisms.
