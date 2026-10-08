#!/usr/bin/env bash
# Full-length hg38 Alu copies from the UCSC RepeatMasker table, as FASTA with the UCSC label in the header.
# usage: extract_alu.sh hg38.fa rmsk.txt.gz OUT.fa [MINLEN=250] [MAXLEN=330]
# rmsk.txt.gz: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz (columns: genoName=6, genoStart=7,
# genoEnd=8, strand=10, repName=11, repClass=12, repFamily=13; check against the table's .sql file).
# Needs bedtools, seqkit. Record the date and the rmsk track's library version (UCSC track description page).
set -euo pipefail
FA=$1; RMSK=$2; OUT=$3; MIN=${4:-250}; MAX=${5:-330}
zcat "$RMSK" | awk -F'\t' -v min="$MIN" -v max="$MAX" '$13 == "Alu" && ($8 - $7) >= min && ($8 - $7) <= max {
    printf "%s\t%d\t%d\t%s|%s|%s:%d-%d\t0\t%s\n", $6, $7, $8, $11, $13, $6, $7, $8, $10 }' > "$OUT.bed"
bedtools getfasta -fi "$FA" -bed "$OUT.bed" -s -name > "$OUT"
echo "$(grep -c '>' "$OUT") copies; labels:"; grep '>' "$OUT" | sed 's/^>//; s/|.*//' | sort | uniq -c | sort -k1,1nr | head -50
