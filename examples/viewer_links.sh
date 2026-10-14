#!/usr/bin/env bash
# Print one ViewAlign / MSA-viewer link per alignment in examples/, pinned to a pushed commit.
#   examples/viewer_links.sh [COMMIT]     (default: HEAD; push it first, otherwise the links do not work)
set -euo pipefail
cd "$(dirname "$0")/.."
C=${1:-$(git rev-parse HEAD)}
git cat-file -e "$C^{commit}" 2>/dev/null || { echo "unknown commit $C" >&2; exit 1; }
C=$(git rev-parse "$C^{commit}")
# the commit must be on the remote, or raw.githubusercontent.com has no such file
git branch -r --contains "$C" 2>/dev/null | grep -q . || echo "warning: $C is not on any fetched remote branch; push it first" >&2
git ls-tree -r --name-only "$C" -- examples/ | grep '\.aln\.fasta$' | while IFS= read -r f; do   # files of that commit, not of the working tree
    t=$(basename "$f" .aln.fasta)
    printf '%s\n  https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/%s/%s&title=%s\n' "$t" "$C" "$f" "$t"
done
