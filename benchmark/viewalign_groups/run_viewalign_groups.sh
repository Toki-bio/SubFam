#!/usr/bin/env bash
# ViewAlign's own group finding (cluster.js: SINEClusterer, "Cluster Now", and characterizeGroups) on SubFam's chunk-consensus alignments.
#   run_viewalign_groups.sh MSAVIEWER_DIR OUTDIR  B2_ALN B2_CHUNKS B2_TRUTH_RMSK B2_TRUTH_V2021  L1TAIL_ALN L1TAIL_CHUNKS L1_TRUTH
# MSAVIEWER_DIR: a checkout of https://github.com/Toki-bio/MSA-viewer (used here: 5d637fb7c2d3b5576599225c8215562293ea371c); node >= 18.
# B2_*: sf50/b2.aln.fasta and b2.chunks.tsv of benchmark/sine_b2/make_b2.sh, b2_truth.tsv, b2_truth_kramerov.tsv.
# L1TAIL_*: SubFam.sh -n 20 -x tail on the last 1000 bp of every L1 copy of benchmark/line_real (l1.fa), l1_truth.tsv as in baseline_vsearch/run_sweeps.sh.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); M=$1; O=$2; mkdir -p "$O"
B2A=$3; B2C=$4; B2R=$5; B2V=$6; LA=$7; LC=$8; LT=$9
{ echo "# B2, 60 chunk consensuses; Cluster Now at the UI defaults, recursion depth 2, columns present in >= 50 % of the rows"
  for mp in 5 2; do
    node "$HERE/cluster_rows.js" "$M" "$B2A" "$O/b2_mp$mp.json" --recurse 2 --min-perfect $mp --min-occupancy 0.5
    echo "## minPerfect $mp"; python3 "$HERE/score_clusters.py" "$O/b2_mp$mp.json" "$B2C" rmsk "$B2R" V2021 "$B2V"
  done
  echo; echo "# B2: exclusive (strict) and enriched (cloudy) columns of a GIVEN partition = chunk majority label (groups below 3 rows are not characterised)"
  python3 - "$B2C" "$B2R" "$B2V" "$O" <<'PY'
import collections, sys
c, r, v, o = sys.argv[1:5]; mem = collections.defaultdict(list)
for l in open(c):
    f = l.split()
    if len(f) >= 2: mem[f[1]].append(f[0])
for nm, tf in (('rmsk', r), ('v2021', v)):
    t = dict(l.split()[:2] for l in open(tf))
    open('%s/b2_rows_%s.tsv' % (o, nm), 'w').write(''.join('%s\t%s\n' % (ch, collections.Counter(t[x] for x in m).most_common(1)[0][0]) for ch, m in sorted(mem.items())))
PY
  for nm in v2021 rmsk; do echo "## groups = chunk majority, $nm"; node "$HERE/cluster_rows.js" "$M" "$B2A" "$O/b2_char_$nm.json" --groups "$O/b2_rows_$nm.tsv"; done
  echo; echo "# L1, last 1000 bp of every copy, 111 chunk consensuses; recursion depth 2, minPerfect 3, columns present in >= 50 % of the rows"
  node "$HERE/cluster_rows.js" "$M" "$LA" "$O/l1tail_occ50.json" --recurse 2 --min-perfect 3 --min-occupancy 0.5
  python3 "$HERE/score_clusters.py" "$O/l1tail_occ50.json" "$LC" rmsk "$LT"
} 2>&1 | tee "$O/viewalign_groups.txt"
