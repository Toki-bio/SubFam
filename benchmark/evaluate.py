#!/usr/bin/env python3
"""Score a set of representative sequences against the simulated truth.

  evaluate.py METHOD REPS.fasta MEMBERS.tsv TRUTH.tsv MASTERS.fasta

MEMBERS.tsv: input id <TAB> group id; each group id must be the name of its
representative in REPS.fasta. Prints one TSV row:

  method  reps  reps_ge10  singleton_frac  purity  median_id  subfam_recovered

reps_ge10        representatives built from >= 10 input sequences
singleton_frac   fraction of input sequences left in groups of size 1
purity           fraction of input sequences whose group majority is their own subfamily
median_id        median identity of the reps_ge10 to their nearest true subfamily consensus
subfam_recovered true subfamilies (of 8) matched at >= 99% identity by a representative
                 whose group mostly comes from that subfamily
"""
import collections
import statistics
import subprocess
import sys
import tempfile


def fasta_ids(path):
    with open(path) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def main():
    method, reps, members, truth, masters = sys.argv[1:6]
    truth = dict(l.split() for l in open(truth))
    groups = collections.defaultdict(list)
    for line in open(members):
        sid, gid = line.split("\t")[:2]
        groups[gid.strip()].append(truth[sid])
    assert sum(map(len, groups.values())) == len(truth), "every input sequence must be in one group"

    majority = {g: collections.Counter(v).most_common(1)[0] for g, v in groups.items()}
    purity = sum(n for _, n in majority.values()) / len(truth)
    singleton = sum(len(v) for v in groups.values() if len(v) == 1) / len(truth)

    rep_ids = fasta_ids(reps)
    assert set(rep_ids) == set(groups), "representative names must match group ids"
    big = [r for r in rep_ids if len(groups[r]) >= 10]

    with tempfile.NamedTemporaryFile("r", suffix=".tsv") as hits:
        subprocess.run(["vsearch", "--usearch_global", reps, "--db", masters, "--id", "0.5",
                        "--iddef", "1", "--strand", "both", "--maxaccepts", "0", "--maxrejects", "0", "--maxhits", "1",
                        "--userout", hits.name, "--userfields", "query+target+id",
                        "--quiet", "--threads", "1"], check=True)
        best = {}
        for line in hits:
            q, t, ident = line.split("\t")
            best[q] = (t, float(ident))

    ids = [best[r][1] if r in best else 0.0 for r in big]
    recovered = {best[r][0] for r in rep_ids
                 if r in best and best[r][1] >= 99.0 and majority[r][0] == best[r][0]}
    print("\t".join([method, str(len(rep_ids)), str(len(big)), f"{singleton:.2f}", f"{purity:.3f}",
                     f"{statistics.median(ids):.1f}" if ids else "NA", f"{len(recovered)}/8"]))


if __name__ == "__main__":
    main()
