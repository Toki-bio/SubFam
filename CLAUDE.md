# Rules for working in this repository

## Alignments are judged by eye, so give clickable links (author's rule, 2026-10-11)
- Whenever a decision, a claim or a figure depends on how an alignment looks, the author must be able to open it and look. Numbers and prose summaries do not replace this.
- Give a viewer link for every alignment mentioned in the reply, as a clickable link in the reply itself (not only a file path), in this form:
  `https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/<COMMIT>/<path>&title=<URL-encoded title>`
  `<COMMIT>` is a pushed commit that contains the file (the commit the reply is about, never a local-only state, never a branch name that may move). Push before linking.
- `examples/viewer_links.sh [COMMIT]` prints the links for all files in `examples/` (default: the current HEAD, which must be pushed).
- Say in the reply which rows to look at (what the author is deciding on), and state that the viewer itself was not opened if it was not.
- When an alignment changes, regenerate and re-send the links; links pinned to older commits can show older row names or older numbers.
- A decision point (keep or drop data, which example, which classification) is presented together with the alignment links that show it.

## Reproducibility and conventions
- Document every run step so that it can be repeated exactly (inputs, seeds, tool versions, commands); each benchmark directory has a README and a script.
- Commit messages end with the two trailer lines the session's attribution reminder gives. Develop and push only on `claude/practical-dijkstra-a0xv8i`; no pull request unless asked.
- No keys or passwords in the repository or in chat. No model identifiers in repository content.
- Short answers (the author is often on a phone); verify from source or text and say what was and was not verified; always give a link when something needed is mentioned; do not flatter.
- Compute runs on the author's servers are not reachable from the cloud container: prepare scripts, do not claim results.
