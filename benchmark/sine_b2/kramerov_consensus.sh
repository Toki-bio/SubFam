#!/usr/bin/env bash
# The mouse B2a-B2e consensus sequences of Vassetzky, Borodulina, Ustyantsev, Kosushkin & Kramerov
# 2021 (Int J Mol Sci 22:9897, doi:10.3390/ijms22189897), taken from Figure S2 of the supplement
# (the only machine-readable place they exist): the supplement zip from Europe PMC, pdftotext,
# then the B2?_Mm rows of the two alignment blocks joined and gap-stripped. IUPAC codes are kept.
#   kramerov_consensus.sh OUT.fa [WORKDIR]
set -euo pipefail
OUT=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); W=${2:-$(mktemp -d)}; mkdir -p "$W"; cd "$W"
[ -s supp.zip ] || curl -sSf -L -o supp.zip 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8466645/supplementaryFiles?includeInlineImage=false'
unzip -o -q supp.zip && unzip -o -q ijms-22-09897-s001.zip
pdftotext -layout "Supplementary Vassetzky/SINE B2. Figs S1-S6. Tables S1-S5/Fig S2.pdf" figS2.txt
python3 -I - figS2.txt "$OUT" <<'EOF'
import re, sys, collections
seqs = collections.OrderedDict()
for l in open(sys.argv[1], encoding='utf-8'):
    m = re.match(r'^(B2[a-e]_Mm)\s+(\S+)\s*$', l)
    if m: seqs.setdefault(m.group(1), []).append(m.group(2))
assert len(seqs) == 5 and all(len(v) == 2 for v in seqs.values()), seqs.keys()
with open(sys.argv[2], 'w') as o:
    for k, v in seqs.items():
        s = ''.join(v).replace('-', ''); assert s.isascii()
        o.write('>KRAMEROV_%s\n%s\n' % (k[:-3], s))
EOF
grep -c '^>' "$OUT"
