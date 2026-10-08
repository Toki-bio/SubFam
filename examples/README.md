# SubFam alignments to open in the viewer

Each link opens the file in [ViewAlign / MSA-viewer](https://toki-bio.github.io/MSA-viewer/) through its `?url=` parameter (the file is fetched from raw.githubusercontent.com, which allows cross-origin requests).
The links point at commit `833b664`, so they keep working while that commit exists; once the files are on the `main` branch, replace the commit id by `main` in the address.
In every file the first rows (`REF_...` or `PRICE_...`) are the true or published consensuses, followed by the SubFam rows in guide-tree order, all aligned together with MAFFT; look for the SubFam row that matches each reference.

| file | what it is | rows | open |
|---|---|---|---|
| sim_simple_middle_n50.aln.fasta | simple simulation, middle family (8 % divergence), 2,000 copies, SubFam `-n 50`; 8 true subfamily masters (`REF_SF1`-`REF_SF8`) + 40 rows | 48 | [open](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/833b6649f6f958fb803dd877e004093703bcff24/examples/sim_simple_middle_n50.aln.fasta&title=SubFam%20simple%20simulation%2C%20middle%20family%2C%20n%3D50) |
| sim_hard_old_n20.aln.fasta | harder simulation (CpG decay, source elements, mixed ages, truncated copies), old family (15 %), SubFam `-n 20`; 8 masters + 100 rows; CpG sites differ from the masters by design | 108 | [open](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/833b6649f6f958fb803dd877e004093703bcff24/examples/sim_hard_old_n20.aln.fasta&title=SubFam%20hard%20simulation%2C%20old%20family%2C%20n%3D20) |
| konkel_alu_n20.aln.fasta | 316 Sanger-sequenced polymorphic Alu loci (GenBank KT305395-KT305737, Konkel et al. 2015) cut to the Alu body, SubFam `-n 20`; Price et al. consensuses AluY, AluYa5, AluYb8, AluYb9 | 19 | [open](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/833b6649f6f958fb803dd877e004093703bcff24/examples/konkel_alu_n20.aln.fasta&title=SubFam%20on%20Konkel%20Alu%20loci%2C%20n%3D20) |
| sim_line_longcopies_n20c.aln.fasta | L1-like simulation (6 kb), rows built from the copies of 2 kb or more with `-n 20 -c`; 8 masters + 15 rows, aligned with `mafft --auto`; rows cover only part of each master | 23 | [open](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/833b6649f6f958fb803dd877e004093703bcff24/examples/sim_line_longcopies_n20c.aln.fasta&title=SubFam%20L1-like%20simulation%2C%20long%20copies) |

How they were made: benchmark/ in this repository (simulate.py, hard/simulate2.py, hard/simulate_line.py, alu_konkel/); the alignments are reproducible from the seeds named in the benchmark READMEs (seed 1 here).
The viewer page itself could not be opened from the environment these links were written in, so the links are untested in the browser; the raw files were checked (HTTP 200, cross-origin allowed).
