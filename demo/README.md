# SubFam demo: the examples of the manuscript, step by step

This folder holds the input copies of every example in the SubFam manuscript (*BMC Bioinformatics*, Software article), so that you can run SubFam on them yourself and get the same rows as in the paper.
Follow the steps in order. Every command is meant to be copied as it is. Each step says what you should see.

| file in `demo/data/` | copies | manuscript | command options |
|---|---|---|---|
| `b2_mm39_chr19_3000.fa` | 3,000 | Example 1, Figure 1 | `-n 50` |
| `b2_mm39_chr19_390.fa` | 390 | Figure 2 | `-n 20` |
| `l1pa_hg38_chr4_2222.fa` | 2,222 | Example 2, Additional file 2 | `-n 20 -c` |
| `alu_konkel2015_316.fa` | 316 | Example 3, Additional file 4 | `-n 20` |
| `b2_mm39_chr1-3_30000.fa` | 30,000 | Example 4, Additional files 5 and 6 | `-n 50` |

Where each file comes from is described at the end of this page.

## Step 1. What you need

- Linux or macOS. On Windows, use WSL (Windows Subsystem for Linux) and run every command inside its Linux terminal.
- conda (Miniforge or Miniconda). Check: `conda --version` prints a version number.
- git. Check: `git --version` prints a version number.
- About 5 GB of free memory for Examples 2 and 4; the other examples need less than 100 MB.

## Step 2. Get SubFam and the demo files

```bash
git clone https://github.com/Toki-bio/SubFam.git
cd SubFam
```

You are now in the folder `SubFam`. All later commands are run from here.

## Step 3. Install the programs SubFam needs

```bash
conda env create -f environment.yml
conda activate subfam
```

This makes a conda environment named `subfam` with MAFFT, gawk, a C compiler and Python. It takes a few minutes, once.
In every new terminal, run `conda activate subfam` again before you use SubFam.

Check that it works:

```bash
./SubFam.sh -v
```

You should see `SubFam 1.5.1`.

**If conda cannot download** (some networks block or throttle its servers; the error is `CondaHTTPError: HTTP 000 CONNECTION FAILED`): on Ubuntu, Debian or WSL install the same programs from the system instead, and skip `conda activate subfam` in the later steps:

```
sudo apt install mafft gawk gcc
```

The results are the same: this was tested with Ubuntu's MAFFT 7.525 and gawk 5.3.2.

## Step 4. Example 1: 3,000 mouse B2 copies

```bash
./SubFam.sh -n 50 -t 8 -o demo_out/ex1 -x b2_3000 demo/data/b2_mm39_chr19_3000.fa
```

`-n 50` puts 50 copies in each chunk, `-t 8` uses 8 processor threads (use the number your computer has), `-o` names the output folder and `-x` the prefix of the output files.

Expected: the count below prints `60` (3,000 copies / 50 = 60 rows), as in the manuscript. It takes about 6 seconds.

Count the rows:

```bash
awk '/^>/{n++} END{print n}' demo_out/ex1/b2_3000.cons.fasta
```

The output folder `demo_out/ex1/` now holds:

| file | what is in it |
|---|---|
| `b2_3000.aln.fasta` | the main result: the chunk consensuses (rows), aligned. Open it in an alignment viewer (Step 9) |
| `b2_3000.cons.fasta` | the same consensuses, not aligned |
| `b2_3000.msf` | the alignment in MSF format |
| `b2_3000.chunks.tsv` | for every input copy: its id, the row it went into, its strand |
| `b2_3000.chunkstats.tsv` | for every row: how many copies, and how well they agree with each other |
| `b2_3000.log` | the SubFam version, the command and the settings |

The rows are named `b2_3000_01`, `b2_3000_02`, ... in the order SubFam found; neighboring numbers are similar.
The copy ids start with their RepeatMasker class (for example `B2_Mm1a_chr19_15922648_15922841`), so you can see in `b2_3000.chunks.tsv` which classes went into which row.

## Step 5. Figure 2: 390 mouse B2 copies

```bash
./SubFam.sh -n 20 -t 8 -o demo_out/fig2 -x b2_390 demo/data/b2_mm39_chr19_390.fa
awk '/^>/{n++} END{print n}' demo_out/fig2/b2_390.cons.fasta
```

Expected: `20`. These 20 rows are panel (c) of Figure 2. The figure itself, with the copies and the rows in one alignment, is drawn by `benchmark/sine_b2/figure_chaos_to_order.py`.

## Step 6. Example 2: 2,222 human L1PA copies, most of them 5′-truncated

```bash
./SubFam.sh -n 20 -c -t 8 -o demo_out/ex2 -x l1pa demo/data/l1pa_hg38_chr4_2222.fa
awk '/^>/{n++} END{print n}' demo_out/ex2/l1pa.cons.fasta
```

`-c` is for truncated copies: a base of a row is called from the copies that reach that position, so the rows reach as far 5′ as their copies do.

Expected: `111`. This is the slow one: about 6 minutes and 4.4 GB of memory with 8 threads, 17 minutes with 2 threads.

## Step 7. Example 3: 316 young human Alu loci

```bash
./SubFam.sh -n 20 -t 8 -o demo_out/ex3 -x alu demo/data/alu_konkel2015_316.fa
awk '/^>/{n++} END{print n}' demo_out/ex3/alu.cons.fasta
```

Expected: `16`. Only 7 of the 16 rows are different sequences: in a young family several chunks give the same consensus (manuscript, Limitations). To count the different rows:

```bash
awk '/^>/{if(s!="")print toupper(s); s=""; next} {s=s $0} END{print toupper(s)}' demo_out/ex3/alu.cons.fasta | sort -u | wc -l
```

This prints `7`.

## Step 8. Example 4: 30,000 mouse B2 copies

```bash
./SubFam.sh -n 50 -t 8 -o demo_out/ex4 -x b2_30000 demo/data/b2_mm39_chr1-3_30000.fa
awk '/^>/{n++} END{print n}' demo_out/ex4/b2_30000.cons.fasta
```

Expected: `600` (588 of them different sequences; use the command of Step 7 on `demo_out/ex4/b2_30000.cons.fasta` to count them). About 4 minutes and 3.6 GB of memory with 8 threads, 12 minutes with 2 threads.

## Step 9. Check that you got exactly the manuscript result

```bash
md5sum demo_out/*/*.cons.fasta
```

(On macOS: `md5 demo_out/*/*.cons.fasta`.) You should see these five checksums:

```
ffabc9786ac57fec4a1da03a6585fa9f  demo_out/ex1/b2_3000.cons.fasta
ef34f2a6031846f20c26e725d77293be  demo_out/ex2/l1pa.cons.fasta
fa66dcef581cfe208ed0626012c11177  demo_out/ex3/alu.cons.fasta
18ef64f9e900f7a2aeec9e5afc9da01b  demo_out/ex4/b2_30000.cons.fasta
bb256d4555bd1191f4e77cbf1adc1860  demo_out/fig2/b2_390.cons.fasta
```

If they match, your rows are the very rows of the manuscript, base for base. The `.aln.fasta` files can differ slightly between MAFFT versions and thread counts (only the 600-row alignment of Example 4 did so in our tests: 263 to 265 columns); the rows themselves do not. Tested with 2 and 8 threads, MAFFT 7.525 and 7.526.

## Step 10. Look at the result

Open any `.aln.fasta` file in an alignment viewer: Jalview, AliView, or ViewAlign in the browser (https://toki-bio.github.io/MSA-viewer/, load the file from your computer).
Read the rows from top to bottom: groups of similar copies appear as blocks of adjacent rows that share the same differences. Lower-case bases are weakly supported: half of the chunk or fewer carry them.

The figures of the manuscript show the same rows grouped by RepeatMasker class and headed by published consensuses. Those grouped alignments are in [`examples/`](../examples/README.md); they are drawn for the reader, and SubFam itself uses no labels.

## Where the input files come from

- `b2_mm39_chr19_3000.fa`: every RepeatMasker B2_Mm1a, B2_Mm1t, B2_Mm2, B3 and B3A copy on chromosome 19 of the mouse genome mm39 (at least 170 bp, no N, on the element strand), 5,972 copies, of which a random 3,000 (seed 1). Made by `benchmark/sine_b2/make_b2.sh`.
- `b2_mm39_chr19_390.fa`: 390 of these copies, a stratified random sample over four groups of Vassetzky et al. 2021 (B2a 120, B2b 30, B2c 120, B2d 120; seed 1; shuffled). Made by `benchmark/sine_b2/chaos_subset.py`; the same file as `benchmark/sine_b2/chaos_subset.fa`.
- `l1pa_hg38_chr4_2222.fa`: every RepeatMasker L1PA2-L1PA8 copy of at least 1 kb on chromosome 4 of the human genome hg38 (no N, on the element strand). Made by `benchmark/line_real/extract_l1.py` from the UCSC `rmsk.txt.gz` table and `chr4.fa.gz`.
- `alu_konkel2015_316.fa`: 316 Sanger-sequenced polymorphic Alu loci of Konkel et al. 2015 (GenBank KT305395-KT305737), cut to the Alu body. The id is `GenBank accession|the closest Price et al. 2004 subfamily`. Made by `benchmark/alu_konkel/prep.py`.
- `b2_mm39_chr1-3_30000.fa`: B2-family copies of mm39 chromosomes 1-3 chosen by the rule of the first file (41,989 copies), a random 30,000 (seed 1). Made by `benchmark/sine_b2/stress_chr123/pipeline.sh`.
