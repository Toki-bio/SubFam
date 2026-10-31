#!/usr/bin/env bash
# Show what "purer under the Vassetzky 2021 labels" means: vsearch 0.60 clusters of the 3,000 B2 copies as rows,
# same format as the SubFam rows (<group>|<V2021 majority>|p<purity>), Vassetzky + Dfam consensuses as heads.
set -euo pipefail
W=/staging/tmp/subfam_b2_purity
mkdir -p "$W"
cd "$W" || exit 1
export TMPDIR=$W/tmp XDG_CACHE_HOME=$W/cache; mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"
T=16
[ -d SubFam ] || git clone -q https://github.com/Toki-bio/SubFam.git
git -C SubFam checkout -q 8f4afe2
TOP=$W/SubFam; B=$TOP/benchmark/sine_b2
mkdir -p dl
# make_b2.sh steps 1-2 (same commands)
[ -s dl/rmsk.txt.gz ] || curl -sSf -o dl/rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz
[ -s dl/chr19.fa.gz ] || curl -sSf -o dl/chr19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr19.fa.gz
if [ ! -s dfam_B2.fa ]; then
    for a in DF000001735 DF000001736 DF000001737 DF000001738 DF000001739; do
        curl -sSf "https://www.dfam.org/api/families/$a/sequence?format=fasta"; echo
    done | grep -v '^$' | sed 's/^>\(DF[0-9.]*\) \(.*\)/>DFAM_\2/' > dfam_B2.fa
fi
python3 -I "$TOP/benchmark/line_real/extract_l1.py" dl/rmsk.txt.gz dl/chr19.fa.gz b2_chr19.fa 170 '^(B2_Mm1a|B2_Mm1t|B2_Mm2|B3|B3A)$'
python3 -I "$B/subsample.py" b2_chr19.fa 3000 1 b2_chr19_3000.fa b2_truth.tsv
# the sample must be the one the committed labels belong to
diff <(grep '^>' b2_chr19_3000.fa | sed 's/^>//' | cut -d' ' -f1 | sort) <(cut -f1 "$B/b2_truth_kramerov.tsv" | sort) > /dev/null && echo "SAMPLE_MATCHES_REPO" || { echo "SAMPLE_DIFFERS"; exit 1; }
cp "$B/b2_truth_kramerov.tsv" "$B/kramerov2021_B2_consensus.fa" .
# vsearch 0.60, exactly as vsearch_sweep.py
vsearch --cluster_fast b2_chr19_3000.fa --id 0.60 --iddef 2 --strand plus --uc vs60.uc --threads $T --quiet
awk -F'\t' '$1=="S"||$1=="H" { split($9, a, " "); print a[1] "\tvs60_" $2 "\t+" }' vs60.uc > vs60_all.tsv
echo "groups: $(cut -f2 vs60_all.tsv | sort -u | wc -l)  singletons: $(cut -f2 vs60_all.tsv | sort | uniq -c | awk '$1==1' | wc -l)"
# clusters of >= 10 copies (the ones scored in Table 2), each turned into a consensus with SubFam itself (one chunk)
awk -F'\t' 'NR==FNR { n[$2]++; next } n[$2] >= 10' vs60_all.tsv vs60_all.tsv > vs60_big.tsv
echo "big groups: $(cut -f2 vs60_big.tsv | sort -u | wc -l)  copies: $(wc -l < vs60_big.tsv)"
rm -f vs60_cons.fa; mkdir -p vs_tmp
for g in $(cut -f2 vs60_big.tsv | sort -u); do
    awk -v g="$g" '$2==g { print $1 }' vs60_big.tsv > vs_tmp/$g.ids
    awk 'NR==FNR { k[$1]; next } /^>/ { split(substr($0, 2), a, " "); p = (a[1] in k) } p' vs_tmp/$g.ids b2_chr19_3000.fa > vs_tmp/$g.fa
    n=$(grep -c '^>' vs_tmp/$g.fa)
    "$TOP/SubFam.sh" -n "$n" -t 4 -o vs_tmp/sf_$g -x c vs_tmp/$g.fa > vs_tmp/$g.log 2>&1
    sed "1s/^>.*/>$g/" vs_tmp/sf_$g/c.cons.fasta >> vs60_cons.fa
done
# rows file in the format build_examples.py takes; group ids sorted by size so the row names read vs60_<cluster>
python3 "$TOP/examples/build_examples.py" b2_vsearch060_V2021.aln.fasta b2_chr19_3000.fa b2_truth_kramerov.tsv \
    vs60_cons.fa vs60_big.tsv --refs kramerov2021_B2_consensus.fa --both --extra-heads dfam_B2.fa --threads $T --tmp build_tmp_v
python3 "$TOP/examples/check_rows.py" b2_vsearch060_V2021.aln.fasta > b2_vsearch060_V2021.rows.tsv
# purity of the big groups under V2021 (Table 2 text says 0.849)
python3 - <<'EOF'
import collections
tr = dict(l.split()[:2] for l in open('b2_truth_kramerov.tsv'))
g = collections.defaultdict(list)
for l in open('vs60_big.tsv'):
    c, k, _ = l.split('\t'); g[k].append(c)
pur = sum(collections.Counter(tr[m] for m in v).most_common(1)[0][1] for v in g.values()) / sum(map(len, g.values()))
print('V2021 purity of groups >= 10: %.3f over %d copies in %d groups' % (pur, sum(map(len, g.values())), len(g)))
EOF
echo ALL_DONE
