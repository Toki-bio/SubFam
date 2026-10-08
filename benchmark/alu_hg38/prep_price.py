#!/usr/bin/env python3
"""Price et al. 2004 supplement (SUBFAMILIES.FINAL.txt) -> FASTA of the 213 consensuses.
usage: prep_price.py SUBFAMILIES.FINAL.txt > price213.fa
Header: >Name|id|copies|P-value (the file is the authors' supplement; keep it outside the repo)."""
import re, sys
n = 0
for l in open(sys.argv[1]):
    m = re.match(r'subfamily\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)(?:\s+\([^)]*\))?\s+([acgtn]+)\s*$', l)
    if m:
        n += 1
        print('>%s|%s|%s|%s\n%s' % (m.group(2), m.group(1), m.group(3), m.group(4), m.group(5).upper()))
sys.stderr.write('%d consensuses\n' % n)
