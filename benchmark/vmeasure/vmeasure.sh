#!/usr/bin/env bash
# Table 2 rescored on ALL copies (singletons are groups of their own): homogeneity, completeness, V-measure,
# next to the old "copies in groups >= 10" and "purity of those" (must reproduce Table 2 = oracle).
set -euo pipefail
W=/staging/tmp/subfam_vmeasure
mkdir -p "$W"
cd "$W" || exit 1
export TMPDIR=$W/tmp XDG_CACHE_HOME=$W/cache; mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" dl groups
T=16
TOP=/staging/tmp/subfam_b2_purity/SubFam; SUB=$TOP/SubFam.sh
B2DIR=/staging/tmp/subfam_b2_purity

# ---- MeShClust2 d7fdb42, built as benchmark/baseline_meshclust/run_meshclust.sh says
MESH=$W/MeShClust2/build/bin/meshclust2
if [ ! -x "$MESH" ]; then
    ( [ -d MeShClust2 ] || git clone -q https://github.com/BioinformaticsToolsmith/MeShClust2.git
      cd MeShClust2 && git checkout -q d7fdb42
      grep -q 'include cstdint' CMakeLists.txt || sed -i 's/set(CMAKE_CXX_FLAGS "/set(CMAKE_CXX_FLAGS "-include cstdint -include limits -include algorithm /' CMakeLists.txt
      mkdir -p build && cd build && cmake .. > ../cmake.log 2>&1 && make -j 8 > ../make.log 2>&1 ) || echo "MESHCLUST_BUILD_FAILED"
    [ -x "$MESH" ] || MESH=$(find "$W/MeShClust2" -type f -name 'meshclust*' -perm -u+x 2>/dev/null | head -1)
fi
echo "meshclust: ${MESH:-none}"

# ---- inputs
# B2: the 3,000 copies of make_b2.sh (rebuilt and checked in subfam_b2_purity)
cp "$B2DIR/b2_chr19_3000.fa" b2.fa; cp "$B2DIR/b2_truth.tsv" b2_rmsk.tsv; cp "$TOP/benchmark/sine_b2/b2_truth_kramerov.tsv" b2_v2021.tsv
# L1: hg38 chr4 L1PA2-8 >= 1 kb (examples/rebuild.sh)
[ -s dl/hg38_rmsk.txt.gz ] || curl -sSf -o dl/hg38_rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
[ -s dl/chr4.fa.gz ] || curl -sSf -o dl/chr4.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr4.fa.gz
python3 -I "$TOP/benchmark/line_real/extract_l1.py" dl/hg38_rmsk.txt.gz dl/chr4.fa.gz l1.fa 1000
grep '^>' l1.fa | sed 's/^>//' | awk -F_ '{ print $0 "\t" $1 }' > l1_rmsk.tsv
# Alu: Konkel et al. 2015 loci from GenBank + Price et al. 2004 supplement (benchmark/alu_konkel/prep.py)
mkdir -p gb
if [ "$(ls gb | wc -l)" -lt 343 ]; then
    for i in $(seq 305395 305737); do
        [ -s gb/KT$i.gb ] || { curl -sSf --retry 5 --retry-delay 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KT$i&rettype=gb&retmode=text" > gb/KT$i.gb; sleep 0.4; }
    done
fi
mkdir -p alu_prep
python3 "$TOP/benchmark/alu_konkel/prep.py" gb SUBFAMILIES.FINAL.txt alu_prep
cp alu_prep/alu_bodies.fa alu.fa
awk -F'|' '/^>/ { print substr($1, 2) "|" $2 "\t" $2 }' alu.fa > alu_price.tsv
echo "copies: b2 $(grep -c '^>' b2.fa)  l1 $(grep -c '^>' l1.fa)  alu $(grep -c '^>' alu.fa)"

# ---- groupings, each as copy<TAB>group
grp_subfam() { "$SUB" "${@:3}" -t $T -o sf_$1 -x $1 $1.fa > sf_$1.log 2>&1; cut -f1,2 sf_$1/$1.chunks.tsv > groups/$1__subfam.tsv; }
grp_subfam b2 x -n 50
grp_subfam l1 x -n 20 -c
grp_subfam alu x -n 20
for s in b2 l1 alu; do grep '^>' $s.fa | sed 's/^>//; s/ .*//' | awk '{ print $1 "\tall" }' > groups/${s}__one.tsv; done
vs() {  # set id
    vsearch --cluster_fast $1.fa --id $2 --iddef 2 --strand plus --uc tmp/v.uc --threads $T --quiet
    awk -F'\t' '$1=="S"||$1=="H" { split($9, a, " "); print a[1] "\t" $2 }' tmp/v.uc > groups/$1__vsearch_$2.tsv
}
ms() {
    [ -n "${MESH:-}" ] && [ -x "$MESH" ] || return 0
    "$MESH" --id $2 --threads $T --output tmp/m.clstr $1.fa > /dev/null 2>&1 || { echo "meshclust $1 $2 failed"; return 0; }
    awk '/^>Cluster/ { c = $2; next } index($0, ">") { x = substr($0, index($0, ">") + 1); sub(/\.\.\..*/, "", x); split(x, d, " "); print d[1] "\t" c }' tmp/m.clstr > groups/$1__meshclust2_$2.tsv
}
for i in 0.60 0.70 0.80 0.90 0.95; do vs b2 $i; ms b2 $i; done
for i in 0.80 0.90 0.95 0.98; do vs l1 $i; ms l1 $i; done
for i in 0.90 0.95 0.97 0.98 0.99 0.995; do vs alu $i; done
for i in 0.90 0.95 0.98 0.99; do ms alu $i; done

# ---- scores
cat > score.py <<'EOF'
import collections, glob, math, sys
def ent(counts):
    n = sum(counts); return -sum(c / n * math.log(c / n) for c in counts if c)
def scores(groups, truth):
    """groups: copy -> group (every copy). Homogeneity/completeness/V (Rosenberg & Hirschberg 2007)."""
    n = len(truth); joint = collections.Counter((truth[c], groups[c]) for c in truth)
    cls = collections.Counter(truth.values()); grp = collections.Counter(groups[c] for c in truth)
    Hc, Hk = ent(cls.values()), ent(grp.values())
    Hc_k = -sum(v / n * math.log(v / grp[k]) for (c, k), v in joint.items())
    Hk_c = -sum(v / n * math.log(v / cls[c]) for (c, k), v in joint.items())
    h = 1 - Hc_k / Hc if Hc else 1.0; cm = 1 - Hk_c / Hk if Hk else 1.0
    v = 2 * h * cm / (h + cm) if h + cm else 0.0
    # the old Table 2 numbers (oracle): share of copies in groups >= 10 and their purity
    members = collections.defaultdict(list)
    for c in truth: members[groups[c]].append(c)
    big = [m for m in members.values() if len(m) >= 10]; cov = sum(map(len, big)) / n
    pur = sum(collections.Counter(truth[c] for c in m).most_common(1)[0][1] for m in big) / max(1, sum(map(len, big)))
    sing = sum(1 for m in members.values() if len(m) == 1)
    return len(members), sing, cov, pur, h, cm, v
sets = [('b2', [('rmsk', 'b2_rmsk.tsv'), ('V2021', 'b2_v2021.tsv')]), ('l1', [('rmsk', 'l1_rmsk.tsv')]), ('alu', [('Price', 'alu_price.tsv')])]
print('set\tlabels\tmethod\tgroups\tsingletons\tcov_ge10\tpurity_ge10\thomogeneity\tcompleteness\tV')
for s, labs in sets:
    for lab, tf in labs:
        truth = dict(l.split()[:2] for l in open(tf))
        def key(f): m = f.split('__')[1][:-4]; return (0 if m == 'one' else 1 if m == 'subfam' else 2, m)
        for f in sorted(glob.glob('groups/%s__*.tsv' % s), key=key):
            g = dict(l.split()[:2] for l in open(f)); miss = [c for c in truth if c not in g]
            for i, c in enumerate(miss): g[c] = '__missing_%d' % i      # a copy a tool dropped = its own group
            r = scores(g, truth)
            print('%s\t%s\t%s%s\t%d\t%d\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f' % (s, lab, key(f)[1], ' (%d missing)' % len(miss) if miss else '', *r))
EOF
python3 score.py | tee vmeasure.tsv
echo ALL_DONE
