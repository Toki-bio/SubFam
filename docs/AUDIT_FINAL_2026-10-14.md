# SubFam pre-publication audit, 2026-10-14

Scope: `SubFam.sh`, `annotate_aligned_by_ssearch_seqkit.sh`, `tests/`, `examples/`, the benchmark drivers and the README.
Start: 1.4.2 (2ff7b62). End: 1.4.4 (7c0713a). Run in a cloud container with MAFFT 7.505, EMBOSS 6.6.0, gawk 5.2.1, mawk, python3 with numpy.

## Verdict

SubFam 1.4.4 is ready for submission. Every defect found was fixed and re-tested. Output is byte-identical to 1.4.2 on the regression input.

## Defects found and fixed

| # | Severity | Where | Defect | Fix, version |
|---|---|---|---|---|
| 1 | fatal | SubFam.sh | 1.4.3 named its checked input `input_checked.fasta` in the work dir; with the prefix `input` the chunk glob matched it and the run stopped ("expected N consensus sequences"). 1.4.2 had the same clash for BOM input only. | input copies in `$WORK/in/`; test with prefix `input`, 1.4.4 |
| 2 | major | SubFam.sh | duplicate ids passed silently: `-e` wrote all copies of an id into one chunk file, `-g` counted one label for all | duplicate id is an error, 1.4.3 |
| 3 | major | SubFam.sh | an empty record entered a chunk (id_min 0, ragged 1) | dropped with a warning, 1.4.3 |
| 4 | minor | SubFam.sh | protein input failed late ("MAFFT failed on chunk"); header without id accepted | clear errors, 1.4.3 |
| 5 | minor | SubFam.sh | `python3 -c` imported a `numpy.py` from the current directory | current directory removed from `sys.path`, 1.4.4 |
| 6 | minor | SubFam.sh | failure message pointed at the work dir the exit trap had deleted | says to rerun with `-K`, 1.4.4 |
| 7 | major | annotate script | strand column always `+` (compared subject, not query coordinates) | `qstart > qend`, 1.4.4 |
| 8 | major | annotate script | no E-value cutoff (an all-N row got a label); `-z 11` shuffles made bit scores differ run to run; no hits against a library of homologs | `-z 3 -E ${ANNOT_EVALUE:-1e-5}`: 3 runs identical, all-N and random rows unlabelled |
| 9 | major | annotate script | exited 1 with no message when nothing matched or the database was missing | explicit checks and messages |
| 10 | minor | annotate script | `python3 -` imported from the current directory; `xargs -n1 -I` warning; `$DB` pasted into `bash -c`; usage showed Clustal input, which seqkit rejects | `-I`, plain `xargs -I`, FASTA in the usage |
| 11 | major | repository | `bench_out/` (an aborted benchmark run with session paths) was committed | removed, in `.gitignore` |
| 12 | major | repository | four documented scripts not executable (`examples/rebuild.sh` gave "Permission denied" from a fresh clone) | mode 755 |
| 13 | minor | stress_chr123 | `pipeline.sh` waited on flag files of this session and read a sibling work dir; `seeds.sh` needed an uncommitted file | both download their own data, take WORKDIR, stop on a failed run |
| 14 | minor | drivers | `run_benchmark.sh` ran 3 min before failing on a missing tool; `multiseed/run.sh` hid SubFam failures; ablation and sensitivity accepted no argument | tool check, `set -e`, `${1:?}` |
| 15 | minor | examples/rebuild.sh | relative `KONKEL`/`L1` silently skipped their example | resolved before `cd`, missing path is an error |
| 16 | cosmetic | docs | README lacked `PREFIX.log`; SCALING.md timed the old full L-INS-i final step; `viewer_links.sh` listed working-tree files | fixed |

## Verified

- `tests/check_equivalence.sh` with `kmer-tree.js`. Consensus against EMBOSS `cons`: 3,600 runs, 0 different. MSF against EMBOSS `seqret`: 80 files, 0 different. k-mer order against ViewAlign `kmer-tree.js`: 20 runs, 0 different.
- `tests/check_inputs.sh` (new): 10 of 10 cases pass.
- Regression: `benchmark/sine_b2/chaos_subset.fa` with `-g` and `-e 1-3` gives byte-identical consensus, alignment, chunk table, chunk statistics, composition and chunk files in 1.4.2 and 1.4.4, under gawk and mawk.
- Option validation: `-n 0`, `-n abc`, `-p 0`, `-p 1.5`, `-k 2`, `-t 0`, an `-e` naming no chunk, a missing or empty file all stop with a clear message. Paths with spaces, gzip, BOM, CRLF, IUPAC codes, lower case, one or two sequences, and `-r -c -a -K`, `-P`, `-I` all run.
- `examples/rebuild.sh` from a fresh clone reproduced `sim_simple_middle_n50` and `sim_hard_old_n20` byte-identical to the committed files.
- The README benchmark table matches `benchmark/results.tsv`. Seeds are fixed in every generator.
- shellcheck at error level is clean on all shell scripts touched; pyflakes shows only unused variables.

## Not verified

- `run_benchmark.sh` end to end: cd-hit-est and mmseqs are not installed here.
- The stress and Alu hg38 drivers were not rerun after the edits; they need about 2 GB of downloads and 4 GB of memory.
- The ViewAlign viewer was checked with Playwright on one example only (it loads; one resource returned 404).

## Note for the bundled copy in SINEderella

The copy of SubFam bundled in SINEderella is 1.2.0. Syncing it to 1.4.4 changes only `input.msf` at SINEderella's defaults, a file SINEderella does not read; `input.clw` is byte-identical at 300 copies. With inputs that are not a multiple of N, 1.4.x cuts round(T/N) equal chunks where 1.2.0 put the remainder in the last chunk, so the files then differ. The SINEderella wrapper has two separate defects (chunks numbered 1000 and up are dropped from `input.clw`; stale chunks of an earlier run are included); they are reported to the author, not fixed here.
