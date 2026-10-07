#!/usr/bin/env python3
"""Simulate a repeat family with known subfamily structure.

A random ancestral sequence evolves along a balanced tree (root -> 2 -> 4 -> 8
subfamilies); every branch fixes a few diagnostic substitutions. Each subfamily
master then spawns copies that decay independently (substitutions + short indels),
like transposon insertions ageing in a genome. Copies are shuffled.

Writes:
  PREFIX.fasta          the simulated copies (input for all methods)
  PREFIX.truth.tsv      copy id -> subfamily
  PREFIX.masters.fasta  true subfamily consensus sequences
"""
import argparse
import random

BASES = "ACGT"


def mutate_sites(seq, k, rng):
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), k):
        seq[pos] = rng.choice([b for b in BASES if b != seq[pos]])
    return "".join(seq)


def decay(seq, rate, rng, indel_share=0.15):
    out = []
    i = 0
    while i < len(seq):
        if rng.random() < rate:
            r = rng.random()
            if r < indel_share / 2:                      # deletion of 1-3 bp
                i += rng.randint(1, 3)
                continue
            if r < indel_share:                          # insertion of 1-3 bp
                out.append("".join(rng.choice(BASES) for _ in range(rng.randint(1, 3))))
                out.append(seq[i])
            else:                                        # substitution
                out.append(rng.choice([b for b in BASES if b != seq[i]]))
        else:
            out.append(seq[i])
        i += 1
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--length", type=int, default=300)
    ap.add_argument("--diag", type=int, default=4, help="diagnostic substitutions per tree branch")
    ap.add_argument("--divergence", type=float, default=0.15, help="per-site decay of each copy from its master")
    ap.add_argument("--sizes", default="600,400,300,250,200,150,60,40", help="copies per subfamily (8 values)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    sizes = [int(x) for x in a.sizes.split(",")]
    assert len(sizes) == 8, "need 8 subfamily sizes"

    level = ["".join(rng.choice(BASES) for _ in range(a.length))]
    for _ in range(3):
        level = [mutate_sites(s, a.diag, rng) for s in level for _ in (0, 1)]
    masters = {f"SF{i + 1}": s for i, s in enumerate(level)}

    copies = []
    for (name, master), n in zip(masters.items(), sizes):
        for j in range(n):
            copies.append((name, decay(master, a.divergence, rng)))
    rng.shuffle(copies)

    # neutral ids, so no method can see the subfamily in the name
    with open(a.prefix + ".fasta", "w") as fa, open(a.prefix + ".truth.tsv", "w") as tr:
        for k, (sf, seq) in enumerate(copies, 1):
            cid = f"copy{k:05d}"
            fa.write(f">{cid}\n{seq}\n")
            tr.write(f"{cid}\t{sf}\n")
    with open(a.prefix + ".masters.fasta", "w") as fa:
        for name, seq in masters.items():
            fa.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    main()
