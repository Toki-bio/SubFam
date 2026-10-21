# Journal options for the SubFam tool paper (2026-10-08; web search results only, no page opened)

Requirement from the owner: free for authors in Uzbekistan; the owner states that the Nature group and BMC journals are free for them. NOT verified here: Uzbekistan's Research4Life group, and whether Springer Nature (BMC, Nature Portfolio) currently waives its charge on that basis
(a 2023 Scholarly Kitchen post says Springer Nature uses World Bank classifications, not Research4Life groups; a university guide says a full BMC waiver in fully open access journals needs all authors in Group A countries). Check before submission:
- Research4Life eligibility lists: https://research4life.org/access/eligibility/
- Research4Life list of open access titles with reduced or no charges, with links to each publisher's waiver policy: https://portal.research4life.org/apclink/getApcLink
- BMC's own waiver page (linked from the portal above) and the journal's current APC.

| journal | why | notes from search results |
|---|---|---|
| **BMC Bioinformatics** (Software article) | scope is novel algorithms and software; the Software template asks for expected wide use and a significant advance shown by direct comparison with related tools (we have VSEARCH, CD-HIT, MMseqs2 comparisons), reviewer access to the tool, free non-commercial availability, source archive and project page; open-source license recommended | https://bmcbioinformatics.biomedcentral.com/submission-guidelines/aims-and-scope ; the software criteria were read from sister BMC journals' pages (same template), not from BMC Bioinformatics' own page |
| **Mobile DNA** (BMC) | the transposable-element journal; two of the Alu papers read in this project were published there (Ahmed 2013, Kryatova 2017); scope includes biocomputational tools (per a CIRAD guide) | a Software article type was not confirmed; aggregator sites list APC about EUR 2,390 / USD 2,890 (not authoritative) |
| Scientific Reports (Nature Portfolio) | broad, but a tool paper must be justified by the benchmark | not checked |

Not Nature/BMC, mentioned only so they are not forgotten: Bioinformatics and NAR Genomics and Bioinformatics (Oxford), GigaScience, Journal of Open Source Software (short software paper; fit not checked).

Before any submission: choose a LICENSE (the README says all rights are reserved until one is added; BMC wants free non-commercial availability and recommends open source), archive a release (Zenodo DOI), check the author list and the name "SubFam" for clashes.

Update: Nikita is not an author; he is thanked in the acknowledgements (README updated). License: MIT (LICENSE added, copyright line Sergei Kosushkin, 2026).

## 2026-10-21

- **SubFam 1.5.1** (6115958): the k-mer ordering is checked against free memory (about 4·n² bytes; SubFam stops with a message
  when it does not fit) instead of the old advice to use `-P` above 20,000 sequences. On real Alu (30,000 copies of hg38 chr1)
  PartTree chunks are less pure than k-mer chunks (0.400 against 0.472, `benchmark/alu_hg38/README.md`), and the k-mer ordering
  ran there in a few GB. Input ids starting with `_R_` are refused (GLM audit 284); the `-a` no-consensus warning also goes to the
  log (GLM audit 362). Tests: `tests/check_inputs.sh` (3 new cases), `test_robust.sh` 26/26, `test_kmer_order.sh` 292/292.
- **Sparse k-mer ordering (`-S`)**, the DeepSeek design tested by AutoClaw on KIT (`/data/W/toki/subfam_sparse_test/`, handoffs in
  `C:\work\SubFam\`): not adopted. It loses to `-P` and to the exact ordering on every real-data test (alu30k 0.378 against 0.399 and
  0.466) and measured only the ordering step, never consensus quality.
- **Manuscript v7** (`docs/MANUSCRIPT_BMC_v7.md`): stale passages corrected (Example 3, Availability, tests, final-alignment modes),
  the memory need of the k-mer ordering and the cost of `-P` stated, abstract 324 words; Figure 2 final
  (`docs/figures/fig2_chaos_to_order`, its MAFFT alignment on one thread so that it is reproducible) and Supplementary Figure 3
  (`benchmark/figure_supp3_simulation.py`); 24 numbered references, each checked against PubMed or the tool page
  (`docs/REFERENCES_BMC_v7.md`). Waiting on the author: prior-use papers, declarations, name spelling, waiver.
