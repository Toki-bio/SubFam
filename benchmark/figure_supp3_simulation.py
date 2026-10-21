#!/usr/bin/env python3
"""figure_supp3_simulation.py OUT_STEM - Supplementary Figure 3 (writes OUT_STEM.png at 300 dpi and OUT_STEM.pdf).

a  Simulated families (multiseed/results.tsv, 8 seeds per age): subfamily sources recovered (of 8) by VSEARCH at each identity
   threshold (mean over seeds, range shaded) against SubFam -n 20 and -n 50 (horizontal lines, mean; seeds that differ in brackets).
b  The number of rows each method returns (median over seeds, log scale).
c  Real B2 copies (sine_b2/sensitivity.tsv): chunk purity against the RepeatMasker class over N (k = 6) and over k (N = 50),
   with four random samples of the copies at the default.
Inputs are read relative to this script; no tool is run."""
import csv, os, statistics, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1]
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "legend.fontsize": 6, "pdf.fonttype": 42})

R = list(csv.DictReader(open(os.path.join(HERE, "multiseed", "results.tsv")), delimiter="\t"))
rec = lambda r: int(r["recovered"].split("/")[0])
ages = [("young", "young (3 %)"), ("middle", "middle (8 %)"), ("old", "old (15 %)")]
ids = sorted({float(r["method"].split()[1]) for r in R if r["method"].startswith("VSEARCH")})

fig = plt.figure(figsize=(6.7, 5.4))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.05], hspace=0.62, wspace=0.32)
col_sf = {"SubFam -n 20": "#1b7837", "SubFam -n 50": "#762a83"}
for j, (age, title) in enumerate(ages):
    rows = [r for r in R if r["scenario"] == age]
    vs = {i: [r for r in rows if r["method"] == "VSEARCH %.2f" % i] for i in ids}
    ax, bx = fig.add_subplot(gs[0, j]), fig.add_subplot(gs[1, j])
    m = [statistics.mean(rec(r) for r in vs[i]) for i in ids]
    ax.fill_between(ids, [min(rec(r) for r in vs[i]) for i in ids], [max(rec(r) for r in vs[i]) for i in ids], color="#bbbbbb", lw=0)
    ax.plot(ids, m, "o-", color="#333333", ms=2.5, lw=1, label="VSEARCH (mean, range)")
    bx.plot(ids, [statistics.median(int(r["reps"]) for r in vs[i]) for i in ids], "o-", color="#333333", ms=2.5, lw=1, label="VSEARCH")
    for meth, c in col_sf.items():
        s = [r for r in rows if r["method"] == meth]
        ax.axhline(statistics.mean(rec(r) for r in s), color=c, lw=1.2, ls="--", label="%s (mean)" % meth)
        bx.axhline(statistics.median(int(r["reps"]) for r in s), color=c, lw=1.2, ls="--", label=meth)
    ax.set_ylim(-0.3, 8.5); ax.set_yticks(range(0, 9, 2)); ax.set_title(title)
    bx.set_yscale("log"); bx.set_ylim(0.7, 3000)
    for x in (ax, bx):
        x.set_xlim(0.63, 0.97); x.set_xticks([0.65, 0.75, 0.85, 0.95]); x.grid(alpha=0.25, lw=0.4)
    bx.set_xlabel("VSEARCH identity threshold")
    if j == 0:
        ax.set_ylabel("sources recovered (of 8)"); bx.set_ylabel("rows returned")
        ax.text(-0.32, 1.12, "a", transform=ax.transAxes, fontweight="bold", fontsize=9)
        bx.text(-0.32, 1.12, "b", transform=bx.transAxes, fontweight="bold", fontsize=9)
    if j == 2:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

S = list(csv.DictReader(open(os.path.join(HERE, "sine_b2", "sensitivity.tsv")), delimiter="\t"))
def pick(n=None, k=None, p="0.36", sample="seed1"):
    for r in S:
        f = dict(x.split("=") for x in r["setting"].split() if "=" in x)
        if r["sample"] == sample and f.get("p") == p and (n is None or f.get("n") == str(n)) and (k is None or f.get("k") == str(k)):
            return float(r["purity_rmsk"])
    raise KeyError((n, k, sample))
cx, kx = fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])
Ns, Ks = [20, 30, 50, 100], [4, 5, 6, 8, 10]
cx.plot(Ns, [pick(n=n, k=6) for n in Ns], "o-", color="#333333", ms=3, lw=1)
seeds = [pick(n=50, k=6, sample="seed%d" % s) for s in (1, 2, 3, 4)]
cx.plot([50] * 4, seeds, "_", color="#c0392b", ms=8, mew=1.2, label="4 random samples, N = 50")
cx.set_xticks(Ns); cx.set_xlabel("N (copies per chunk), k = 6")
cx.set_ylabel("chunk purity"); cx.legend(frameon=False, loc="lower left")
kx.plot(Ks, [pick(n=50, k=k) for k in Ks], "o-", color="#333333", ms=3, lw=1)
kx.set_xticks(Ks); kx.set_xlabel("k (k-mer size), N = 50")
for x in (cx, kx):
    x.set_ylim(0.70, 0.84); x.grid(alpha=0.25, lw=0.4)
cx.text(-0.32, 1.22, "c", transform=cx.transAxes, fontweight="bold", fontsize=9)
cx.set_title("3,000 real B2 copies (mm39 chr19)", loc="left")
fig.savefig(out + ".png", dpi=300, bbox_inches="tight"); fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png", out + ".pdf", file=sys.stderr)
