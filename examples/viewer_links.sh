#!/usr/bin/env bash
# Print one ViewAlign / MSA-viewer link per alignment in examples/, pinned to a pushed commit.
#   examples/viewer_links.sh [COMMIT]     (default: HEAD; push it first, otherwise the links do not work)
set -euo pipefail
cd "$(dirname "$0")/.."
C=${1:-$(git rev-parse HEAD)}
git cat-file -e "$C^{commit}" 2>/dev/null || { echo "unknown commit $C" >&2; exit 1; }
for f in examples/*.aln.fasta; do
    t=$(basename "$f" .aln.fasta)
    printf '%s\n  https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/SubFam/%s/%s&title=%s\n' "$t" "$C" "$f" "$t"
done
