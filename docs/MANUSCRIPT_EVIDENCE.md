# SubFam manuscript: what the evidence supports (working document, 2026-10-08)

No draft exists. This lists each statement a manuscript would want to make, what backs it today, and what is missing.
"Supported" means a result exists that can be reproduced from this repository; the command or file is named.
Updated as results arrive; entries marked PENDING are running or queued.

| # | Statement | Status | Evidence / what is missing |
| --- | --- | --- | --- |
| 1 | SubFam reduces a large set of related copies to an alignment of chunk consensuses with no threshold other than the chunk size | Supported (design) | `SubFam.sh`; the chunk consensus is a plurality call, equal to EMBOSS `cons` on 3,600 random alignments (`tests/check_equivalence.sh`) |
| 2 | The k-mer ordering is byte-identical across implementations and fast enough for 30,000+ copies | Supported | `kmer_order.c` vs the Python/numpy version: 292 synthetic cases and real SINE sets (`tests/test_kmer_order.sh`); 30,000 copies in 150 s on 16 threads, 1,218 s on 1 thread; identical to ViewAlign's `kmer-tree.js` on 10 comparisons |
| 3 | Whole-pipeline run time and memory against input size, per ordering method | PENDING | `benchmark/scaling/run_scaling.sh` (needs a quiet machine) |
| 4 | The output is reproducible: same input, any thread count, any directory, any day | Supported | `tests/test_robust.sh` (24 cases, determinism over 1/8/8/16 threads on 3,000 real copies, byte-identical); MSF header carries no path and honours `SOURCE_DATE_EPOCH` |
| 5 | SubFam handles awkward input without silent loss | Supported | `tests/test_robust.sh`: single record, fewer records than a chunk, duplicate ids, spaces and pipes in ids, lower case / N / IUPAC, empty record, CRLF, wrapped lines, identical sequences, gapped input, sequences shorter than k, mixed strands, protein (refused with a message) |
| 6 | Chunks recover the subfamilies of a simulated family with a known tree | PENDING | `tools/sim` in sinederella (`coseg-compare` branch), running on KIT: 13 scenarios, purity and ARI per chunk size |
| 7 | Chunks recover the coarse Y lineages of real Alu loci | Supported, limited | 316 Konkel loci, chunk purity 0.918 (n 20) against labels derived from similarity to Price consensuses (`benchmark/alu_konkel`); the labels are not independent of the method |
| 8 | SubFam chunks are not a partition into subfamilies; they over-split by design | Supported | Four hand-curated SINE sets: best ARI over all copies 0.19-0.27 at n 50, falling as n shrinks (sinederella `docs/COMPARISON_LOG.md`) |
| 9 | The SINEderella route (chunks + peel) matches COSEG | Mixed, not supported in general | ARI over all copies, COSEG vs chunks + peel: saq 0.447 vs 0.195, ccr 0.316 vs 0.087, teu 0.429 vs 0.430, dmo 0.424 vs 0.416 (500-900 copies per set, too few for the peel); the 16,000-copy simulation is the fair test and is PENDING |
| 10 | SubFam works for LINEs and truncated copies (`-c`) | Not supported | only a simulated check of the coverage-relative plurality; no LINE data tested |
| 11 | "Subfamily identification" | Not supported | a chunk is not a subfamily; the call is made by a person or the peel |

## Reviewer questions to prepare for

* Why is COSEG (Price et al. 2004) not the reference method, and what was the comparison? (Entry 9 and the comparison log.)
* What independent truth was used? (Only the simulation has true labels; hand-curated groups are one expert's call; the Alu labels are similarity-derived.)
* Do chunks smooth minority diagnostic changes? Yes: a majority-rule consensus over 50 copies cannot carry a change present in a minority of each chunk (SUBFAMILY_METHOD.md section 7 of SINE_discriminator).
* How does the result depend on the chunk size? (Entry 8; the simulation sweeps n.)
