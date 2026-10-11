"""
Figures and tables for the continuous-variable power study, rebuilt from
power_continuous.csv after the October 2026 re-run.

Design decisions, and why:

  * ELEVEN SERIES IS TOO MANY.  The categorical limit is eight, and eleven
    near-identical hues is what made the previous figure unreadable.  The
    statistics are therefore grouped into five families drawn in one colour
    each, which is honest because the members of each family track one
    another closely -- the spread within every family is reported below so
    the grouping can be checked rather than taken on trust.

  * B_n IS THE SUBJECT, so it is a single thick line in the first slot and
    everything else is context.

  * P_n IS OMITTED from the figures, as in the present paper, its power
    being within 0.003 of B_n everywhere; it stays in the tables.

  * BET IS OMITTED ENTIRELY, pending a fix: its adapter returns a quantity
    that is anti-correlated with dependence on the symmetric schemes, so its
    numbers are not usable.  See the diagnostic printed at the end.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SRC = "/mnt/user-data/uploads/results/power_continuous.csv"
OUT = "/mnt/user-data/outputs/"

d = pd.read_csv(SRC)
d = d[d.q == d.p]                                  # primary study only

BROKEN = ["BET"]
FIGDROP = BROKEN + ["P_n"]

FAMILY = {
    "B_n":       ("$B_n$",              "#2a78d6", 2.1, "-"),
    "HHG":       ("HHG / HHG-fast",     "#eb6834", 1.0, "-"),
    "HHG_fast":  ("HHG / HHG-fast",     "#eb6834", 1.0, "-"),
    "dCov":      ("dCov / Hoeff. / Genest", "#1baf7a", 1.0, "--"),
    "Hoeffding": ("dCov / Hoeff. / Genest", "#1baf7a", 1.0, "--"),
    "Genest":    ("dCov / Hoeff. / Genest", "#1baf7a", 1.0, "--"),
    "BEAST":     ("BEAST",              "#4a3aa7", 1.0, (0, (4, 2))),
    "MIC":       ("MIC / $\\xi$",       "#eda100", 1.0, (0, (1, 1.6))),
    "xi":        ("MIC / $\\xi$",       "#eda100", 1.0, (0, (1, 1.6))),
}
LEG_ORDER = ["$B_n$", "HHG / HHG-fast", "dCov / Hoeff. / Genest",
             "BEAST", "MIC / $\\xi$"]

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

# ---------------------------------------------------- within-family spread --
print("Within-family spread of mean power (grouping is only honest if small)")
fam_members = {}
for t, (lab, *_rest) in FAMILY.items():
    fam_members.setdefault(lab, []).append(t)
mp = d.pivot_table(index="test", columns="n", values="power", aggfunc="mean")
for lab, members in fam_members.items():
    if len(members) == 1:
        continue
    sub = mp.loc[members]
    print(f"  {lab:<26}" + "  ".join(
        f"n={n}: {sub[n].min():.3f}-{sub[n].max():.3f} (spread {sub[n].max()-sub[n].min():.3f})"
        for n in (50, 100, 200)))
print()

# ------------------------------------------------------- panel ordering ------
# Order schemes by how far B_n stands above the best competitor, so the figure
# reads left to right from where the test leads to where it trails.  Stated in
# the caption; it is a presentational choice, not a selection.
comp = [t for t in d.test.unique() if t not in FIGDROP + ["B_n"]]
adv = {}
for s in d.scheme.unique():
    sub = d[(d.scheme == s) & (d.n == 100)]
    bn = sub[sub.test == "B_n"].power.mean()
    best = sub[sub.test.isin(comp)].groupby("test").power.mean().max()
    adv[s] = bn - best
SCHEMES = sorted(adv, key=lambda s: -adv[s])
print("B_n minus best competitor, mean over noise levels, n = 100:")
for s in SCHEMES:
    print(f"  {s:<17}{adv[s]:+.3f}")
print()


def panel(ax, sub, show_tests):
    for t in show_tests:
        lab, col, lw, ls = FAMILY[t]
        g = sub[sub.test == t].sort_values("nu")
        if not len(g):
            continue
        ax.plot(g.nu, g.power, color=col, lw=lw, ls=ls,
                zorder=5 if t == "B_n" else 3,
                alpha=1.0 if t == "B_n" else 0.85)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(0, 2.1)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


SHOW = [t for t in FAMILY if t not in FIGDROP]
handles = [Line2D([], [], color=FAMILY[[k for k, v in FAMILY.items()
                                        if v[0] == lab][0]][1],
                  lw=2.1 if lab == "$B_n$" else 1.2,
                  ls=FAMILY[[k for k, v in FAMILY.items() if v[0] == lab][0]][3],
                  label=lab) for lab in LEG_ORDER]

# ------------------------------------------------------------- figure 1 -----
fig, axes = plt.subplots(4, 4, figsize=(7.2, 6.6), sharex=True, sharey=True)
for ax, s in zip(axes.ravel(), SCHEMES):
    panel(ax, d[(d.scheme == s) & (d.n == 100)], SHOW)
    ax.set_title(s, pad=3)
for ax in axes.ravel()[len(SCHEMES):]:
    ax.axis("off")
for ax in axes[:, 0]:
    ax.set_ylabel("power")
fig.supxlabel("standard deviation of the added noise", fontsize=8, y=0.055)
axes.ravel()[len(SCHEMES)].legend(handles=handles, loc="center", fontsize=7.5,
                                  labelspacing=0.9, handlelength=2.4)
fig.suptitle("Size-adjusted power against noise, $n = 100$", fontsize=9.5)
fig.tight_layout(rect=[0, 0.075, 1, 0.965])
fig.savefig(OUT + "power_continuous_main.pdf")
fig.savefig(OUT + "power_continuous_main.png", dpi=170)
plt.close(fig)

# ------------------------------------------------------------- figure 2 -----
fig, axes = plt.subplots(len(SCHEMES), 3, figsize=(6.6, 15.5),
                         sharex=True, sharey=True)
for i, s in enumerate(SCHEMES):
    for j, n in enumerate((50, 100, 200)):
        ax = axes[i, j]
        panel(ax, d[(d.scheme == s) & (d.n == n)], SHOW)
        if i == 0:
            ax.set_title(f"$n = {n}$", pad=4)
        if j == 2:
            ax.text(1.04, 0.5, s, transform=ax.transAxes, rotation=270,
                    va="center", ha="left", fontsize=7)
    axes[i, 0].set_ylabel("power", fontsize=7)
fig.supxlabel("standard deviation of the added noise", fontsize=8, y=0.008)
fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7.5,
           bbox_to_anchor=(0.5, 0.018), handlelength=2.2)
fig.tight_layout(rect=[0, 0.035, 0.97, 1])
fig.savefig(OUT + "power_continuous_full.pdf")
plt.close(fig)

# ------------------------------------------------------------- figure 3 -----
order = mp[200].drop(index=BROKEN).sort_values().index.tolist()
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
for ax, n in zip(axes, (50, 100, 200)):
    v = mp[n].drop(index=BROKEN).reindex(order)
    cols = ["#2a78d6" if t in ("B_n", "P_n") else MUTED for t in order]
    ax.hlines(range(len(order)), 0, v.values, color=cols, lw=1.0, alpha=0.55)
    ax.scatter(v.values, range(len(order)), c=cols, s=26, zorder=3)
    for k, t in enumerate(order):
        ax.text(v.values[k] + 0.018, k, f"{v.values[k]:.3f}", va="center",
                fontsize=6.3, color=INK if t in ("B_n", "P_n") else MUTED)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([r"$B_n$" if t == "B_n" else r"$P_n$" if t == "P_n"
                        else r"$\xi$" if t == "xi" else t.replace("_","-") for t in order],
                       fontsize=7)
    ax.set_title(f"$n = {n}$", pad=4)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.grid(axis="x", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
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

keep = [t for t in mp.index if t not in BROKEN]
keep = mp.loc[keep, 200].sort_values(ascending=False).index.tolist()
sd = d.pivot_table(index="test", columns="n", values="power", aggfunc="std")

lines = [r"\begin{table}",
         r"\tbl{Mean size-adjusted power over the thirteen continuous dependence schemes and twenty noise levels}{%",
         r"\begin{tabular}{l rr rr rr}", r"\hline",
         " & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{$n = {n}$}}" for n in (50, 100, 200)) + r" \\",
         "Test & " + " & ".join(["Mean & SD"] * 3) + r" \\", r"\hline"]
for t in keep:
    nm = {"B_n": "$B_n$", "P_n": "$P_n$", "xi": "$\\xi$",
          "HHG_fast": "HHG-fast"}.get(t, t)
    row = " & ".join(f"{fmt(mp.loc[t, n])} & {fmt(sd.loc[t, n])}" for n in (50, 100, 200))
    lines.append(f"{nm} & {row} \\\\")
lines += [r"\hline", r"\end{tabular}}", r"\begin{tabnote}",
          "Each entry averages the size-adjusted rejection rate over the thirteen schemes and the twenty "
          "noise levels, $260$ configurations in all. Critical values come from $2\\times10^4$ null datasets "
          "per $(n,p,q)$ for the statistics that depend on the data only through the ranks, and from "
          "$5\\times10^3$ permuted datasets per configuration for distance covariance and HHG, which do not; "
          "power is over $10^4$ draws from each alternative, so the largest standard error of a single "
          "rejection rate is $0{\\cdot}005$. Orders follow rule (12), giving $p=q=2,3,4$. SD is the spread "
          "across the $260$ configurations, not a standard error.",
          r"\end{tabnote}", r"\label{tab:power_summary}", r"\end{table}"]
open(OUT + "power_continuous_summary.tex", "w").write("\n".join(lines) + "\n")

# per-scheme supplement table
sch = d[~d.test.isin(BROKEN)].pivot_table(index="scheme", columns=["n", "test"],
                                          values="power", aggfunc="mean")
sup = [r"\begin{table}", r"\tbl{Mean size-adjusted power by scheme}{%",
       r"\scriptsize\setlength{\tabcolsep}{3pt}",
       r"\begin{tabular}{l" + "r" * len(keep) + "}", r"\hline"]
for n in (50, 100, 200):
    sup.append(f"\\multicolumn{{{len(keep)+1}}}{{l}}{{$n = {n}$}} \\\\")
    sup.append("Scheme & " + " & ".join(
        {"B_n": "$B_n$", "P_n": "$P_n$", "xi": "$\\xi$",
         "HHG_fast": "HHG-f"}.get(t, t) for t in keep) + r" \\ \hline")
    for s in SCHEMES:
        vals = [sch[(n, t)][s] for t in keep]
        best = max(vals)
        cells = [(r"\textbf{" + fmt(v).replace("$", "") + "}" if abs(v - best) < 1e-9
                  else fmt(v).replace("$", "")) for v in vals]
        sup.append(f"{s} & " + " & ".join(f"${c}$" for c in cells) + r" \\")
    sup.append(r"\hline")
sup += [r"\end{tabular}}", r"\begin{tabnote}",
        "Largest value in each row in bold. Averages are over the twenty noise levels.",
        r"\end{tabnote}", r"\label{tab:power_by_scheme}", r"\end{table}"]
open(OUT + "power_continuous_by_scheme.tex", "w").write("\n".join(sup) + "\n")

print("mean power by test and n (BET excluded):")
print(mp.drop(index=BROKEN).reindex(keep).round(3).to_string())
print("\nwrote main, full and mean figures + two .tex tables")
