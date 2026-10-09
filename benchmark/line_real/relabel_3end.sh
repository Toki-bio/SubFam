#!/usr/bin/env bash
# Re-label the extracted L1 copies by their last 900 bp, best hit to the Dfam 4.0 L1PA2-8 3'-end consensuses
# (dfam4_L1PA2-8_3end.fa, fetched from https://www.dfam.org/api/families/DF0000003{39..45}/sequence?format=fasta), and tabulate against rmsk.
# usage: relabel_3end.sh l1.fa OUTDIR      -> OUTDIR/hits3p.tsv  (the cross-tabulation is printed by the python block in README.md section 2)
set -euo pipefail; HERE=$(cd "$(dirname "$0")" && pwd); mkdir -p "$2"
python3 -I - "$1" "$2/l1_3p.fa" <<'P'
import sys
n=None;s=[]
out=open(sys.argv[2],'w')
for l in open(sys.argv[1]):
    l=l.strip()
    if l.startswith('>'):
        if n: out.write('>%s\n%s\n'%(n,''.join(s)[-900:]))
        n=l[1:];s=[]
    else: s.append(l)
out.write('>%s\n%s\n'%(n,''.join(s)[-900:]))
P
vsearch --usearch_global "$2/l1_3p.fa" --db "$HERE/dfam4_L1PA2-8_3end.fa" --id 0.5 --iddef 2 --strand plus --maxaccepts 0 --maxrejects 0 --top_hits_only \
  --userout "$2/hits3p.tsv" --userfields query+target+id+qcov+tcov --threads 4 --quiet
