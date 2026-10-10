# Line-by-line review of SubFam.sh v1.2.0

Review pass, 2026-10-10. Line numbers refer to SubFam.sh at commit 682301a. **[verified]** = reproduced
by running the script on crafted input; the two must-fix bugs 1.1 and 1.2 were reproduced a second time
independently (mawk 1.3.4: `runaway regular expression`; `100*0.07 = 7.000000000000001`,
`20*0.36 = 7.199999999999999`, `60*0.36 = 21.599999999999998`).

## What it does well
- The weighted-Jaccard matrix is **exact** (brute-force check on 44,850 pairs incl. repetitive, <k and
  empty sequences: 0 mismatches) **[verified]**, computed with a "sum of min = Σ_t [a≥t][b≥t]" dense-BLAS
  trick and no per-pair Python loop.
- UPGMA is genuine average linkage (size-weighted update, line 250) with a lazy nearest-neighbour cache,
  O(n²) in practice (5,000×300 bp: distances 7 s, UPGMA 3 s) **[verified]**.
- Edge cases behave: 1 sequence, 2 sequences, N > total, sequences < k, empty record, CRLF, duplicate /
  empty / `|` / `>` / space IDs all complete with exit 0 and sane output **[verified]**. Output is byte-
  identical between `-t 1` and `-t 4` and between runs **[verified]**.
- Option validation is good; `set -euo pipefail` plus xargs exit 123 means a failing chunk aborts the run;
  the consensus awk is faithful to EMBOSS `cons` and has an equivalence test.

## 1. Correctness bugs

**1.1 mawk cannot parse the chunk splitter (README says mawk is supported) — crash.** Line 440:
`c = int(i++ / n); if (c >= nc) c = nc - 1` — mawk reads `/ n); ...` after `++` as a regex and exits 2
with no outputs **[verified]**. Fix: `c = int(i / n); i++`.

**1.2 Plurality threshold off by one for many (n, p) pairs — floating-point ceil.** Line 457 and 347:
`x = n * p; (x == int(x)) ? x : int(x) + 1`. `100*0.07` → 8 instead of 7; same for `-p 0.14 -n 50/100`,
`-p 0.28 -n 25/50/100`, `-p 0.36 -n 20/60` **[verified]**. Default 0.36 at n=50 is exact (18.0), but the
`-c` path (`need`) hits whatever `ncov` a column has, so under `-c` some columns get the wrong threshold.
`examples/build_examples.py:42` (`need = -(-frac * ncov // 1)`) has the same defect. Fix:
`int(n*p + 1 - 1e-9)` or integer arithmetic.

**1.3 Prefix derived before basename → hidden/garbled output names.** Line 112:
`PREFIX=${PREFIX:-$(basename "${IN%.*}")}`. `./d/noext` → PREFIX `""` → hidden dot-files **[verified]**;
`data.v2/seqs` → `data`; `x.fa.gz` → `x.fa`. Fix: `b=$(basename "$IN"); PREFIX=${PREFIX:-${b%.*}}`.

**1.4 Prefix with glob metacharacters silently produces empty results, exit 0.** Lines 453/462 use
`find -name "${NAME}_*.fasta"`; `-x 'a[1]'` matches nothing, outputs are 0 bytes, "Job completed"
**[verified]**. `-x a/b` dies in `mktemp`. Fix: bash glob or `printf '%s\0'`, validate `-x`, die if no
`.cons` file was produced.

**1.5 gzip input is rejected by accident.** Lines 419–420 `grep -c '^>'` on the gzip binary gave 0 →
"no FASTA records" **[verified]**; another payload could pass and then `open()` raises
`UnicodeDecodeError`. Support `.gz` (`zcat` into `$WORK`) or detect and die clearly.

**1.6 numpy version not checked.** Line 105 tests `import numpy` only; line 155 needs
`sliding_window_view` (numpy ≥ 1.20); line 206 `intersect1d(return_indices=True)` (≥ 1.15).

**1.7 `THREADS=$(nproc)` under `set -e`.** Without `nproc` (macOS) the script dies before `-h`.
Use `nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`.

**1.8 `-a` inconsistency:** no-consensus positions are `n` (lowercase, via `setcase`), low-coverage
positions under `-c` are `N` (line 349 continues before the case rule) **[verified]**. README says `N`.

**1.9 IDs beginning with `_R_`** are reported as reverse-complemented (lines 443–444). Minor.

Non-issues checked: duplicate IDs, IUPAC/lowercase/N, chunk of 1 or 2, N > TOTAL, very short sequences.

## 2. Algorithmic issues

**Weighted Jaccard (163–200).** True multiset (Ružička) Jaccard on k-mer counts, exact, float32. k-mers
are taken over A/C/G/T/U only, so an `N` is *deleted* and its neighbours joined into a spurious k-mer
(line 152); masked copies look slightly more similar than they are. Canonical k-mers (`-r`) correct.
Cost: the count-level loop (181) builds a dense `n × block` matrix for *every* level t up to the maximum
k-mer count; poly-A tails (count ≈ 55) made the 5,000-sequence matrix 5× slower (7.4 → 37 s)
**[verified]**. Compacting the column index per level removes this.

**UPGMA (212–264).** Genuine UPGMA, global minimum merge. Ties break to the lowest row index, then lowest
column: deterministic but input-order dependent. Leaf order: lower-index cluster first, each side flipped
so the closest *end leaves* meet (241–246) — a greedy "closest ends" heuristic, not optimal leaf ordering
(Bar-Joseph); consistent with kmer-tree.js and pinned by the equivalence test. Describe it as such.

**Chunking (421, 437–447).** Consecutive blocks of the leaf order, never clades. `NCHUNK = floor(TOTAL/N)`,
remainder joins the *last* chunk, so it holds N…2N−1 sequences (120 seqs, `-n 50` → 50 + 70)
**[verified]**. Prefer `NCHUNK = round(TOTAL/N)` or spread the remainder.

**Consensus (339–367).** EMBOSS semantics: winner = residue of the first sequence with the maximum
EDNAFULL column score; a gap scores 0 (360) and beats any residue with a negative score — 2G/2T/1gap
gives no call although G reaches plurality 2 **[verified]**. Gap-majority columns yield `N`, stripped by
default, so the consensus is gap-free; with `-a`, deletion columns become `n`, indistinguishable from
"no consensus". Plurality counts EDNAFULL-positive matches (`mt`), so IUPAC `R` counts toward A and G.
`setcase = n/2` lowercases any base supported by ≤ 50 % of sequences — with p = 0.36 a large fraction of
the output, undocumented. **`-c`** (334–349, 359) matches its description. `build_examples.py:30–45`
picks by raw count, not EDNAFULL score, so reference rows can differ from SubFam's own rule at columns
with IUPAC codes or a gap veto. **0.36** is justified nowhere in README/usage; `docs/SCALING.md:218`
shows it is 18/50 inherited from SINEderella.

**Options [verified from code]:** `-x` prefix for files, temp dir and names (spaces → `_` in names).
`-r`: canonical k-mers for the tree **and** greedy orientation toward the predecessor; the frame is the
first sequence of the order, so the output can come out reverse-complemented relative to the input
majority **[verified]**. `-m`: `mafft --retree 0 --reorder`. `-P`: adds `--parttree`. `-a`: skip the N
strip. `-K`: skip the trap. `--ep 0.123` (469) is MAFFT's default offset — a no-op.

## 3. Scalability
- Memory: dense float32 n×n (178): 400 MB at 10⁴, 40 GB at 10⁵, impossible at 10⁶. Sparse table (169–176)
  ≈ 4 int64 arrays of total k-mer length: 10⁵ × 1 kb ≈ 3–4 GB transient.
- Time: 5,000 seqs 7 s; ~10⁴ ≈ 30 s; 10⁵ ≈ hours; UPGMA minutes at 2×10⁴. The Python is single-threaded
  except BLAS, which ignores `-t` (set `OPENBLAS_NUM_THREADS=$THREADS`).
- Final alignment uses L-INS-i (469): impractical beyond a few thousand consensuses. Fall back to `--auto`
  when `NCHUNK` is large.
- Temp dir removed on *any* exit (116): a failed run leaves nothing without `-K`.
- The split awk never `close()`s chunk files (446); non-gawk awks with low `ulimit -n` break.

## 4. Robustness / UX
- `-p 0.5x` is accepted (awk numeric prefix) **[verified]**.
- The `awk | awk` pipeline (459) has no `pipefail` in `sh`, so a failing consensus awk yields an empty
  `.cons` silently; a failed MAFFT surfaces only as xargs exit 123 without the chunk name.
- Dependency check misses numpy version and mafft version (≥ 7 for `--adjustdirection`/`--thread`).
- Reproducibility: no RNG; MAFFT deterministic here across thread counts **[verified at 120 seqs]**; add
  `--threadit 0` to 469 if byte-reproducibility matters.
- No run log of parameters/version in the output directory.
- `aln.fasta` rows are in MAFFT `--reorder` order (469), not chunk order; `cons.fasta` and `chunks.tsv`
  are in chunk order.

## 5. README ↔ code mismatches
- "python3 with numpy" — needs numpy ≥ 1.20.
- "tested with gawk and mawk" — mawk fails (1.1).
- "Leftover sequences (fewer than N) join the last chunk" — true, but the last chunk is up to 2N−1.
- "A base is called when at least `-p` × chunk size of the sequences agree" — actually EDNAFULL-positive
  matches of the top-scoring residue, gaps score 0 and can veto, IUPAC codes count as partial matches.
- `-a` "keep … as `N`" — `n` except `-c` low-coverage positions.
- Undocumented: lowercase = supported by ≤ 50 %; names `PREFIX_01…` zero-padded; `-r` frame; final
  alignment row order is MAFFT's.
- `docs/SCALING.md:166` mentions `-c MIN`; the code has no MIN argument (`MINCOV=3` hard-coded).
- Usage "best at all divergences tested" and 0.36 have no justification in README.

## 6. Prioritised fixes

**Must-fix before publication**
1. Line 440: `c = int(i / n); i++`.
2. Lines 457 and 347 (and `build_examples.py:42`): ceil with epsilon / integer arithmetic.
3. Line 112: basename first, then strip the extension; refuse an empty prefix.
4. Lines 453/462: no `find -name` with an unescaped prefix; validate `-x`; die if no `.cons` files.
5. Line 105: check `numpy >= 1.20` and document it.

**Should-fix**
6. Line 421: balance chunks; document the size range.
7. Lines 349/364–365: one symbol for both no-call kinds under `-a`, or document.
8. Line 181: compact k-mer columns per count level; export `OPENBLAS_NUM_THREADS=$THREADS`.
9. Line 469: choose MAFFT mode by `NCHUNK`; add `--threadit 0`; drop `--ep 0.123`.
10. Line 459: run the two awks in sequence or check `.cons` is non-empty; print the chunk name on failure.
11. Line 24: portable CPU count; lines 419–420: detect gzip.
12. Line 446: `close()` the previous chunk file.

**Nice**
13. Orient `-r` output toward the majority strand.
14. Write a run log into `$OUTDIR`; per-chunk progress.
15. Document `setcase`, `--reorder` order, the gap veto, the origin of 0.36; make `MINCOV` an option.
16. Align `build_examples.py:plurality()` with the awk rule, or state that reference rows use a count rule.
