#!/usr/bin/env python3
"""Harder repeat-family simulation than benchmark/simulate.py.

Adds to the simple model (balanced tree, 8 subfamilies, independent decay):
  * CpG-rich ancestor and CpG hypermutability: a CpG site mutates (C>T, G>A) at CPG x the rate of other sites.
  * Source elements: each subfamily has 1-5 active sources (master + a few private mutations) with skewed
    activity, instead of every copy coming straight from the subfamily master.
  * Age mixture: each copy has its own divergence, drawn around the scenario mean (gamma, CV 0.35);
    optionally a different mean per subfamily (--ages), so one dataset mixes young and old subfamilies.
  * Transitions are twice as likely as transversions; indels are 10 % of mutation events.
  * Truncation: a fraction of the copies are 5'-truncated fragments (min length + exponential).
Writes PREFIX.fasta, PREFIX.truth.tsv (copy -> subfamily), PREFIX.masters.fasta (subfamily consensus before source-private mutations).
"""
import argparse
import random

BASES = "ACGT"
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def make_ancestor(length, cpg_sites, rng):
    s = [rng.choice(BASES) for _ in range(length)]
    for pos in rng.sample(range(length - 1), cpg_sites):
        s[pos], s[pos + 1] = "C", "G"
    return s


def substitute(base, rng):
    if rng.random() < 2 / 3:                      # transitions twice as likely as each pair of transversions
        return TRANSITION[base]
    return rng.choice([b for b in BASES if b != base and b != TRANSITION[base]])


def mutate_sites(seq, k, rng):
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), k):
        seq[pos] = substitute(seq[pos], rng)
    return "".join(seq)


def cpg_mask(seq):
    m = [False] * len(seq)
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            m[i] = m[i + 1] = True
    return m


def decay(seq, rate, cpg_mult, rng, indel_share=0.10):
    mask = cpg_mask(seq)
    out = []
    i = 0
    while i < len(seq):
        b = seq[i]
        if mask[i] and cpg_mult > 1:
            if rng.random() < min(0.9, rate * cpg_mult):
                out.append("T" if b == "C" else "A")      # CpG decay: C>T, G>A
                i += 1
                continue
        if rng.random() < rate:
            r = rng.random()
            if r < indel_share / 2:
                i += rng.randint(1, 3)
                continue
            if r < indel_share:
                out.append("".join(rng.choice(BASES) for _ in range(rng.randint(1, 3))))
                out.append(b)
            else:
                out.append(substitute(b, rng))
        else:
            out.append(b)
        i += 1
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--cpg-sites", type=int, default=24, help="CpG dinucleotides in the ancestor")
    ap.add_argument("--cpg", type=float, default=6.0, help="CpG mutation rate multiplier (1 = none)")
    ap.add_argument("--diag", type=int, default=4, help="diagnostic substitutions per tree branch")
    ap.add_argument("--divergence", type=float, default=0.08, help="mean per-site decay of a copy from its source (non-CpG sites)")
    ap.add_argument("--ages", default="", help="8 comma-separated per-subfamily mean divergences (overrides --divergence)")
    ap.add_argument("--sizes", default="600,400,300,250,200,150,60,40")
    ap.add_argument("--sources", type=int, default=5, help="maximum active sources per subfamily")
    ap.add_argument("--private", type=float, default=2.0, help="mean private mutations of a source relative to its subfamily master")
    ap.add_argument("--trunc-frac", type=float, default=0.30, help="fraction of copies that are 5'-truncated")
    ap.add_argument("--trunc-min", type=int, default=100)
    ap.add_argument("--trunc-mean", type=int, default=90, help="mean extra length of a truncated copy beyond trunc-min")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    sizes = [int(x) for x in a.sizes.split(",")]
    ages = [float(x) for x in a.ages.split(",")] if a.ages else [a.divergence] * 8
    assert len(sizes) == 8 and len(ages) == 8

    level = ["".join(make_ancestor(a.length, a.cpg_sites, rng))]
    for _ in range(3):
        level = [mutate_sites(s, a.diag, rng) for s in level for _ in (0, 1)]
    masters = {"SF%d" % (i + 1): s for i, s in enumerate(level)}

    copies = []
    for (name, master), n, age in zip(masters.items(), sizes, ages):
        k = rng.randint(1, a.sources)
        sources = []
        for r in range(k):
            npriv = 0 if r == 0 else max(1, int(rng.expovariate(1 / a.private)) + 1)   # the dominant source is the master
            sources.append(mutate_sites(master, npriv, rng) if npriv else master)
        weights = [(r + 1) ** -1.5 for r in range(k)]
        for _ in range(n):
            src = rng.choices(sources, weights)[0]
            d = rng.gammavariate(1 / 0.35 ** 2, age * 0.35 ** 2)       # mean age, CV 0.35
            c = decay(src, d, a.cpg, rng)
            if rng.random() < a.trunc_frac:
                keep = min(len(c), a.trunc_min + int(rng.expovariate(1 / a.trunc_mean)))
                c = c[len(c) - keep:]
            copies.append((name, c))
    rng.shuffle(copies)
    with open(a.prefix + ".fasta", "w") as fa, open(a.prefix + ".truth.tsv", "w") as tr:
        for k, (sf, seq) in enumerate(copies, 1):
            fa.write(">copy%05d\n%s\n" % (k, seq))
            tr.write("copy%05d\t%s\n" % (k, sf))
    with open(a.prefix + ".masters.fasta", "w") as fa:
        for name, seq in masters.items():
            fa.write(">%s\n%s\n" % (name, seq))


if __name__ == "__main__":
    main()
