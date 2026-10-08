#!/usr/bin/env python3
"""L1-like simulation: long, mostly 5'-truncated copies, some with a 5'-inverted part and internal deletions.

Reuses the CpG-aware decay of simulate2.py. A 6 kb master with a CpG-rich 5' region (first 900 bp) and a CpG-poor remainder evolves along a balanced tree
(8 subfamilies). Each branch fixes DIAG substitutions, half of them in the 5' region (visible only in long copies), half in the last 1.5 kb (visible in most copies).
Copies: source elements as in simulate2.py, per-copy age (gamma, CV 0.35), CpG decay, then
  * 5' truncation: length = trunc-min + exponential(trunc-mean), capped at the full length; FULL of the copies stay full length
  * a fraction INV of the copies have the 5' part (a random 30-90 % of the copy) reverse-complemented (twin priming)
  * a fraction DEL have one internal deletion of 50-500 bp
Writes PREFIX.fasta, PREFIX.truth.tsv, PREFIX.masters.fasta, PREFIX.copyinfo.tsv (copy, subfamily, length, full-length?, inverted?)."""
import argparse, importlib.util, os, random

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("simulate2", os.path.join(here, "simulate2.py"))
s2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s2)
BASES = "ACGT"


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_master(length, rng, utr=900, cpg_utr=70, cpg_rest=60):
    s = [rng.choice(BASES) for _ in range(length)]
    for lo, hi, k in ((0, utr, cpg_utr), (utr, length, cpg_rest)):
        for pos in rng.sample(range(lo, hi - 1), k):
            s[pos], s[pos + 1] = "C", "G"
    return "".join(s)


def mutate_region(seq, k, lo, hi, rng):
    seq = list(seq)
    for pos in rng.sample(range(lo, hi), k):
        seq[pos] = s2.substitute(seq[pos], rng)
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--length", type=int, default=6000)
    ap.add_argument("--diag", type=int, default=20, help="diagnostic substitutions per branch (half 5' region, half 3' 1.5 kb)")
    ap.add_argument("--divergence", type=float, default=0.10)
    ap.add_argument("--cpg", type=float, default=6.0)
    ap.add_argument("--sizes", default="600,400,300,250,200,150,60,40")
    ap.add_argument("--sources", type=int, default=5)
    ap.add_argument("--private", type=float, default=3.0)
    ap.add_argument("--full", type=float, default=0.03, help="fraction of copies kept full length")
    ap.add_argument("--trunc-min", type=int, default=200)
    ap.add_argument("--trunc-mean", type=int, default=900)
    ap.add_argument("--inv", type=float, default=0.15)
    ap.add_argument("--dele", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    sizes = [int(x) for x in a.sizes.split(",")]
    L = a.length

    level = [make_master(L, rng)]
    for _ in range(3):
        nxt = []
        for s in level:
            for _ in (0, 1):
                s = mutate_region(s, a.diag // 2, 0, 900, rng) if False else s
                c = mutate_region(s, a.diag // 2, 0, 900, rng)
                c = mutate_region(c, a.diag - a.diag // 2, L - 1500, L, rng)
                nxt.append(c)
        level = nxt
    masters = {"SF%d" % (i + 1): s for i, s in enumerate(level)}

    copies = []; info = []
    for (name, master), n in zip(masters.items(), sizes):
        k = rng.randint(1, a.sources); sources = []
        for r in range(k):
            sources.append(master if r == 0 else s2.mutate_sites(master, max(1, int(rng.expovariate(1 / a.private)) + 1), rng))
        weights = [(r + 1) ** -1.5 for r in range(k)]
        for _ in range(n):
            src = rng.choices(sources, weights)[0]
            d = rng.gammavariate(1 / 0.35 ** 2, a.divergence * 0.35 ** 2)
            c = s2.decay(src, d, a.cpg, rng)
            full = rng.random() < a.full
            if not full:
                keep = min(len(c), a.trunc_min + int(rng.expovariate(1 / a.trunc_mean)))
                c = c[len(c) - keep:]
            if rng.random() < a.dele and len(c) > 800:
                dl = rng.randint(50, 500); st = rng.randint(100, len(c) - dl - 100)
                c = c[:st] + c[st + dl:]
            inv = rng.random() < a.inv and len(c) > 400
            if inv:
                b = int(len(c) * rng.uniform(0.3, 0.9)); c = revcomp(c[:b]) + c[b:]
            copies.append((name, c, full, inv))
    rng.shuffle(copies)
    with open(a.prefix + ".fasta", "w") as fa, open(a.prefix + ".truth.tsv", "w") as tr, open(a.prefix + ".copyinfo.tsv", "w") as ci:
        for k, (sf, seq, full, inv) in enumerate(copies, 1):
            cid = "copy%05d" % k
            fa.write(">%s\n%s\n" % (cid, seq)); tr.write("%s\t%s\n" % (cid, sf))
            ci.write("%s\t%s\t%d\t%d\t%d\n" % (cid, sf, len(seq), full, inv))
    with open(a.prefix + ".masters.fasta", "w") as fa:
        for name, seq in masters.items():
            fa.write(">%s\n%s\n" % (name, seq))


if __name__ == "__main__":
    main()
