"""
Figures and tables for the continuous-variable power study, rebuilt from
power_continuous.csv after the October 2026 top-up that added the order
scan and repaired BET.

What changed since the previous build, and why the design changed with it:

  * B_scan IS NOW THE SUBJECT.  It has the highest mean power at every
    sample size, and its margin over B_n grows with n (+0.005, +0.043,
    +0.075), which is the reverse of B_n's own pattern.  So the scan gets
    the thick line and B_n becomes its reference, drawn in the same hue
    and a lighter weight: they are the same statistic, scanned and fixed.

  * BET IS BACK.  The signed-statistic bug is fixed, so it rejoins the
    figures and tables.  It is grouped with BEAST, both being binary
    expansion methods from the same authors, which keeps the hue count
    at five where the categorical limit is eight.

  * P_n STAYS OUT OF THE FIGURES, within 0.003 of B_n everywhere, and
    stays in the tables.

  * TWELVE SERIES IS STILL TOO MANY FOR TWELVE COLOURS.  The statistics
    are grouped into families drawn in one hue each; the within-family
    spread is printed below so the grouping can be checked rather than
    believed.  Palette validated for colour-vision deficiency separation:
    worst adjacent pair deltaE 9.2 deutan, 27.6 normal, all five inside
    the lightness band.  Each family also carries its own dash pattern,
    so identity never rests on hue alone.

21 October 2026: the two data-driven rank tests of Kallenberg and Ledwina
(1999), T_S2 (diagonal Legendre components, d(n) = 4) and V (mixed
components, d(n) = 2), joined the study.  They take the next slot of the
categorical order, magenta #e87ba4.  They do NOT track one another -- T_S2
collapses on the heteroscedastic, quadratic and cubic schemes where V does
not -- so they share the hue but not the line style, and each has its own
legend entry.  Palette re-validated with the six hues: adjacent pairs pass
(worst CVD deltaE 9.2 deutan, normal 19.6); all-pairs, magenta-orange is the
weakest pair (normal deltaE 12.9), as yellow-orange already was (13.7), and
both are separated by dash pattern.  Paths are now relative to this file:
it reads power_continuous.csv beside it and writes into ./figs.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import os
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "power_continuous.csv")
OUT = os.path.join(HERE, "figs", "")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(SRC)
d = d[d.q == d.p]                                   # primary study only
NS = (50, 100, 200)

FIGDROP = ["P_n"]                                   # tables only

#   name                     legend label            hue        lw   dash
FAMILY = {
    "B_scan":    ("Order scan",                "#2a78d6", 2.2, "-"),
    "B_n":       ("$B_n$, rule (12)",          "#2a78d6", 1.0, (0, (5, 2))),
    "HHG":       ("HHG / HHG-fast",            "#eb6834", 1.0, "-"),
    "HHG_fast":  ("HHG / HHG-fast",            "#eb6834", 1.0, "-"),
    "dCov":      ("dCov / Hoeff. / Genest",    "#1baf7a", 1.0, (0, (4, 1.6))),
    "Hoeffding": ("dCov / Hoeff. / Genest",    "#1baf7a", 1.0, (0, (4, 1.6))),
    "Genest":    ("dCov / Hoeff. / Genest",    "#1baf7a", 1.0, (0, (4, 1.6))),
    "BET":       ("BET / BEAST",               "#4a3aa7", 1.0, (0, (3, 1.4, 1, 1.4))),
    "BEAST":     ("BET / BEAST",               "#4a3aa7", 1.0, (0, (3, 1.4, 1, 1.4))),
    "MIC":       ("MIC / $\\xi$",              "#eda100", 1.0, (0, (1, 1.6))),
    "xi":        ("MIC / $\\xi$",              "#eda100", 1.0, (0, (1, 1.6))),
    "KL_V":      ("KL $V$",                    "#e87ba4", 1.0, (0, (7, 1.6))),
    "KL_TS2":    ("KL $T_{S2}$",               "#e87ba4", 1.0, (0, (1.6, 1.3, 1.6, 1.3, 5, 1.3))),
}
LEG_ORDER = ["Order scan", "$B_n$, rule (12)", "KL $V$", "KL $T_{S2}$",
             "HHG / HHG-fast", "dCov / Hoeff. / Genest", "BET / BEAST",
             "MIC / $\\xi$"]

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dcdbd7"
plt.rcParams.update({
    "font.family": "serif", "font.size": 8,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.titlesize": 8, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

PRETTY = {"B_n": "$B_n$", "P_n": "$P_n$", "xi": "$\\xi$",
          "HHG_fast": "HHG-fast", "B_scan": "Order scan",
          "KL_V": "KL $V$", "KL_TS2": "KL $T_{S2}$"}
mp = d.pivot_table(index="test", columns="n", values="power", aggfunc="mean")
sd = d.pivot_table(index="test", columns="n", values="power", aggfunc="std")

# ---------------------------------------------------- within-family spread --
print("Within-family spread of mean power (grouping is honest only if small)")
fam = {}
for t, (lab, *_r) in FAMILY.items():
    fam.setdefault(lab, []).append(t)
for lab, members in fam.items():
    if len(members) == 1:
        continue
    sub = mp.loc[members]
    print(f"  {lab:<26}" + "  ".join(
        f"n={n}: {sub[n].min():.3f}-{sub[n].max():.3f}" for n in NS))

# ------------------------------------------------------- panel ordering -----
# Ordered by how far the scan stands above the best competitor at n = 100, so
# the figure reads left to right from where it leads to where it trails.  A
# presentational choice, stated in the caption; nothing is selected on it.
comp = [t for t in d.test.unique() if t not in ("B_scan", "B_n", "P_n")]
adv = {}
for s in d.scheme.unique():
    sub = d[(d.scheme == s) & (d.n == 100)]
    adv[s] = (sub[sub.test == "B_scan"].power.mean()
              - sub[sub.test.isin(comp)].groupby("test").power.mean().max())
SCHEMES = sorted(adv, key=lambda s: -adv[s])
print("\nscan minus best competitor, mean over noise levels, n = 100:")
for s in SCHEMES:
    print(f"  {s:<17}{adv[s]:+.3f}")


def panel(ax, sub, tests):
    for t in tests:                       # B_scan last so it sits on top
        lab, col, lw, ls = FAMILY[t]
        g = sub[sub.test == t].sort_values("nu")
        if len(g):
            ax.plot(g.nu, g.power, color=col, lw=lw, ls=ls,
                    zorder=6 if t == "B_scan" else 5 if t == "B_n" else 3,
                    alpha=1.0 if t in ("B_scan", "B_n") else 0.85)
    ax.set_ylim(-0.03, 1.03); ax.set_xlim(0, 2.1)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


SHOW = [t for t in FAMILY if t not in FIGDROP and t != "B_scan"] + ["B_scan"]
handles = []
for lab in LEG_ORDER:
    k = next(k for k, v in FAMILY.items() if v[0] == lab)
    handles.append(Line2D([], [], color=FAMILY[k][1],
                          lw=2.2 if lab.endswith("scan") else 1.2,
                          ls=FAMILY[k][3], label=lab))

# ------------------------------------------------------------ figure 1 ------
fig, axes = plt.subplots(4, 4, figsize=(7.2, 7.5), sharex=True, sharey=True)
for ax, s in zip(axes.ravel(), SCHEMES):
    panel(ax, d[(d.scheme == s) & (d.n == 100)], SHOW)
    ax.set_title(s, pad=3)
## Thirteen panels in a 4 x 4 grid leaves the last row almost empty, and
## sharex then strips the tick labels from every column whose bottom panel
## is not in that row -- three of the four, with no numbers on the axis at
## all.  Put them back on the lowest panel actually drawn in each column.
flat = axes.ravel()
for col in range(4):
    lowest = max(i for i in range(len(SCHEMES)) if i % 4 == col)
    flat[lowest].tick_params(labelbottom=True)
for ax in flat[len(SCHEMES):]:
    ax.axis("off")
for ax in axes[:, 0]:
    ax.set_ylabel("power")
fig.supxlabel("standard deviation of the added noise", fontsize=8, y=0.028)
flat[len(SCHEMES)].legend(handles=handles, loc="center", fontsize=7,
                          labelspacing=0.85, handlelength=2.6)
fig.suptitle("Size-adjusted power against noise, $n = 100$", fontsize=9.5)
fig.tight_layout(rect=[0, 0.045, 1, 0.968])
fig.savefig(OUT + "power_continuous_main.pdf")
fig.savefig(OUT + "power_continuous_main.png", dpi=170)
plt.close(fig)

# ------------------------------------------------------------ figure 2 ------
## Thirteen rows by three columns is 15.5 inches tall, which cannot be put
## on an A4 page at \textwidth: LaTeX would scale it to 6.5 inches wide and
## 15.3 tall and run it off the bottom.  Split it in two, seven rows then
## six, each of which fits a page with its caption.
for part, block in (("a", SCHEMES[:7]), ("b", SCHEMES[7:])):
    h = 1.15 * len(block) + 0.95
    fig, axes = plt.subplots(len(block), 3, figsize=(6.5, h),
                             sharex=True, sharey=True)
    for i, s in enumerate(block):
        for j, n in enumerate(NS):
            ax = axes[i, j]
            panel(ax, d[(d.scheme == s) & (d.n == n)], SHOW)
            if i == 0:
                ax.set_title(f"$n = {n}$", pad=4)
            if j == 2:
                ax.text(1.04, 0.5, s, transform=ax.transAxes, rotation=270,
                        va="center", ha="left", fontsize=7)
        axes[i, 0].set_ylabel("power", fontsize=7)
    ## The axis label sits under the axes and the legend under that; at the
    ## previous settings the two were written on top of one another.
    yl = 0.66 / h
    fig.supxlabel("standard deviation of the added noise", fontsize=8, y=yl)
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.2,
               bbox_to_anchor=(0.5, 0.004), handlelength=2.6,
               columnspacing=1.4)
    fig.tight_layout(rect=[0, 1.00 / h, 0.97, 1])
    fig.savefig(OUT + f"power_continuous_full_{part}.pdf")
    fig.savefig(OUT + f"power_continuous_full_{part}.png", dpi=150)
    plt.close(fig)

# ------------------------------------------------------------ figure 3 ------
order = mp[200].sort_values().index.tolist()
fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
for ax, n in zip(axes, NS):
    v = mp[n].reindex(order)
    cols = ["#2a78d6" if t in ("B_scan", "B_n", "P_n") else MUTED for t in order]
    ax.hlines(range(len(order)), 0, v.values, color=cols, lw=1.0, alpha=0.55)
    ax.scatter(v.values, range(len(order)),
               c=cols, s=[34 if t == "B_scan" else 24 for t in order], zorder=3)
    for k, t in enumerate(order):
        ax.text(v.values[k] + 0.018, k, f"{v.values[k]:.3f}", va="center",
                fontsize=6.3,
                color=INK if t in ("B_scan", "B_n", "P_n") else MUTED)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([PRETTY.get(t, t) for t in order], fontsize=7)
    ax.set_title(f"$n = {n}$", pad=4)
    ax.set_xlim(0, 1.0); ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.grid(axis="x", color=GRID, lw=0.4); ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
axes[1].set_xlabel("mean size-adjusted power over the thirteen schemes")
fig.tight_layout()
fig.savefig(OUT + "power_continuous_mean.pdf")
fig.savefig(OUT + "power_continuous_mean.png", dpi=170)
plt.close(fig)


# --------------------------------------------------------------- tables -----
def fmt(v):
    return f"$0{{\\cdot}}{int(round(v*1000)):03d}$" if v < 1 else "$1{\\cdot}000$"


keep = mp[200].sort_values(ascending=False).index.tolist()
lines = [r"\begin{table}",
         r"\tbl{Mean size-adjusted power over the thirteen continuous dependence "
         r"schemes and twenty noise levels}{%",
         r"\begin{tabular}{l rr rr rr}", r"\hline",
         " & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{$n = {n}$}}" for n in NS) + r" \\",
         "Test & " + " & ".join(["Mean & SD"] * 3) + r" \\", r"\hline"]
for t in keep:
    row = " & ".join(f"{fmt(mp.loc[t, n])} & {fmt(sd.loc[t, n])}" for n in NS)
    lines.append(f"{PRETTY.get(t, t)} & {row} \\\\")
lines += [
    r"\hline", r"\end{tabular}}", r"\begin{tabnote}",
    "Each entry averages the size-adjusted rejection rate over the thirteen schemes "
    "and the twenty noise levels, $260$ configurations in all. The order scan is the "
    "order scan of \\S\\ref{sec:scan}, which reports the smallest exact $p$-value "
    "among the orders $k = 1,\\ldots,K(n)$ with $K = 4, 5, 6$, calibrated on the "
    "joint null; $B_n$ uses the single order of rule~(12). Critical values come from "
    "$2\\times10^4$ null datasets per $(n,p,q)$ for the statistics that depend on the "
    "data only through the ranks, and from $5\\times10^3$ permuted datasets per "
    "configuration for distance covariance and HHG, which do not; power is over $10^4$ "
    "draws from each alternative, so the largest standard error of a single rejection "
    "rate is $0{\\cdot}005$. Orders follow rule~(12), giving $p=q=2,3,4$, and the binary "
    "expansion depth for BET and BEAST is $d=2,3,3$, matched so that neither family is "
    "favoured by its tuning. SD is the spread across the $260$ configurations, not a "
    "standard error.",
    r"\end{tabnote}", r"\label{tab:power_summary}", r"\end{table}"]
open(OUT + "power_continuous_summary_new.tex", "w").write("\n".join(lines) + "\n")

sch = d.pivot_table(index="scheme", columns=["n", "test"], values="power",
                    aggfunc="mean")
SHORT = {"B_n": "$B_n$", "P_n": "$P_n$", "xi": "$\\xi$", "HHG_fast": "HHG-f",
         "B_scan": "scan", "Hoeffding": "Hoeff", "BEAST": "BEAST"}
sup = [r"\begin{table}", r"\tbl{Mean size-adjusted power by scheme}{%",
       r"\scriptsize\setlength{\tabcolsep}{2.5pt}",
       r"\begin{tabular}{l" + "r" * len(keep) + "}", r"\hline"]
for n in NS:
    sup.append(f"\\multicolumn{{{len(keep)+1}}}{{l}}{{$n = {n}$}} \\\\")
    sup.append("Scheme & " + " & ".join(SHORT.get(t, t) for t in keep) + r" \\ \hline")
    for s in SCHEMES:
        vals = [sch[(n, t)][s] for t in keep]
        best = max(vals)
        sup.append(f"{s} & " + " & ".join(
            ("$\\mathbf{" + fmt(v).replace("$", "") + "}$" if abs(v - best) < 1e-9
             else fmt(v)) for v in vals) + r" \\")
    sup.append(r"\hline")
sup += [r"\end{tabular}}", r"\begin{tabnote}",
        "Largest value in each row in bold; averages are over the twenty noise levels. "
        "Columns are ordered by mean power at $n = 200$, rows by the margin of the scan "
        "over its best competitor at $n = 100$.",
        r"\end{tabnote}", r"\label{tab:power_by_scheme}", r"\end{table}"]
open(OUT + "power_continuous_by_scheme.tex", "w").write("\n".join(sup) + "\n")

# ------------------------------------------- the facts the prose needs ------
print("\n--- mean power by test and n ---")
print(mp.reindex(keep).round(3).to_string())
print("\n--- scan versus rule (12) ---")
for n in NS:
    print(f"  n={n}: scan {mp.loc['B_scan', n]:.3f}  B_n {mp.loc['B_n', n]:.3f}"
          f"  gap {mp.loc['B_scan', n] - mp.loc['B_n', n]:+.3f}")
print("\n--- per-scheme: scan, B_n, best competitor ---")
for s in SCHEMES:
    r = {t: sch[(100, t)][s] for t in keep}
    bc = max((v, t) for t, v in r.items() if t not in ("B_scan", "B_n", "P_n"))
    print(f"  {s:<17} scan {r['B_scan']:.3f}  B_n {r['B_n']:.3f}  "
          f"best {bc[0]:.3f} ({bc[1]})")
print("\n--- how often each competitor beats the scan, over the 13 schemes ---")
for t in keep:
    if t in ("B_scan", "P_n"):
        continue
    w = sum(sch[(n, t)][s] > sch[(n, "B_scan")][s] for s in SCHEMES for n in NS)
    print(f"  {t:<11}{w:>3} of 39 scheme-by-size combinations")
print("\n--- the two Kallenberg-Ledwina tests do not track one another ---")
print((mp.loc["KL_V"] - mp.loc["KL_TS2"]).round(3).to_string())
print("\nwrote three figures and two .tex tables")
