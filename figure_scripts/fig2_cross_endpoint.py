# -*- coding: utf-8 -*-
"""Figure 2 (final numbering; file fig3_cross_endpoint.py).

Fig 3 - cross-endpoint differences D = q95(A) - q95(B), whole-cochlea atlas.

Layout note (BMC compliance): single-column figure, plotted at 170 mm wide so
that nominal font size == effective font size; all text >= 7 pt; height inside
the 225 mm cap. The two panels are stacked vertically.

Panel A: forest plot of the 12 most significant pairwise tests (by MC p),
         with BH-FDR q values and the prior-correction survival flag.
Panel B: lollipop plot of all 84 pairwise tests sorted by MC p (most
         significant at top); the six BH-FDR < 0.05 tests are highlighted.

Input (read-only):
  project/analysis/scdrs_main/whole_cochlea_atlas/compare_pairwise_tests.csv
"""
import csv, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import C_SIG, C_GREY, C_DARK, save, add_panel_label

ATLAS = ("C:/Users/docta/WorkBuddy/纯生信/project/analysis/scdrs_main/"
         "whole_cochlea_atlas")

with open(os.path.join(ATLAS, "compare_pairwise_tests.csv"), newline="") as fh:
    P = list(csv.DictReader(fh))

for r in P:
    r["diff"] = float(r["diff"])
    r["mc_p"] = float(r["mc_p"])
    r["q"] = float(r["qval_BH"])
    r["n"] = int(float(r["n_cell"]))
P.sort(key=lambda r: r["mc_p"])
for i, r in enumerate(P):
    r["rank"] = i + 1
    r["survives"] = r["q"] * 2 < 0.05

sig = [r for r in P if r["q"] < 0.05]

FS_TICK = 7.0
FS_ANNOT = 7.0
FS_LABEL = 8.0
FS_TITLE = 9.0
FS_LEGEND = 7.0

fig = plt.figure(figsize=(6.69, 8.0))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 1.42], hspace=0.40,
                       left=0.205, right=0.985, top=0.945, bottom=0.068)
LIM = (-3.4, 2.6)


def q_text_x(D):
    """Place the q annotation on the emptier side of the point."""
    return D + 0.15 if D < 0.6 else D - 0.15


# ------------------------------ Panel A -----------------------------------
N_TOP = 12
top = P[:N_TOP]
axA = fig.add_subplot(gs[0, 0])
for i, r in enumerate(top):
    col = C_SIG if r["survives"] else (C_DARK if r["q"] < 0.05 else C_GREY)
    axA.plot([LIM[0], r["diff"]], [i, i], color="#D5D5D5", lw=1.0, zorder=1)
    axA.scatter(r["diff"], i, s=58 if r["survives"] else 42, marker="o",
                facecolor=col, edgecolor=C_DARK if r["survives"] else col,
                linewidths=0.7, zorder=3)
    if r["survives"]:
        txt = "q = %.4g\nq\u2032 = %.4g (survives)" % (r["q"], r["q"] * 2)
        axA.annotate(txt, xy=(q_text_x(r["diff"]), i), xycoords="data",
                     fontsize=FS_ANNOT, color=C_DARK, va="center", ha="left",
                     fontweight="bold", linespacing=1.35)
    else:
        axA.annotate("q = %.4g" % r["q"], xy=(q_text_x(r["diff"]), i),
                     xycoords="data", fontsize=FS_ANNOT,
                     color="#444444" if r["q"] >= 0.05 else C_DARK,
                     va="center", ha="left" if r["diff"] < 0.6 else "right")
axA.axvline(0, color=C_DARK, lw=0.9, zorder=2)
axA.set_yticks(np.arange(N_TOP))
axA.set_yticklabels(["%s (n=%s)\n%s vs %s" % (r["cell_type"], f"{r['n']:,}",
                                              r["ep_A"], r["ep_B"])
                     for r in top], fontsize=FS_TICK)
axA.invert_yaxis()
axA.set_xlim(*LIM)
axA.set_ylim(N_TOP - 0.45, -0.85)
axA.set_xlabel("D = q95(A) $-$ q95(B)   (negative: second endpoint $>$ first)",
               fontsize=FS_LABEL)
axA.grid(axis="x", ls=":", lw=0.4, color="#EDEDED", zorder=0)
axA.tick_params(axis="y", length=0)
axA.tick_params(axis="x", labelsize=FS_TICK)
axA.legend(handles=[
    Line2D([], [], marker="o", ls="none", mfc=C_SIG, mec=C_DARK, ms=7,
           label="BH-FDR < 0.05 and q\u2032 = 2q < 0.05"),
    Line2D([], [], marker="o", ls="none", mfc=C_DARK, mec=C_DARK, ms=6,
           label="BH-FDR < 0.05"),
    Line2D([], [], marker="o", ls="none", mfc=C_GREY, mec=C_GREY, ms=6,
           label="BH-FDR \u2265 0.05")],
    loc="center right", bbox_to_anchor=(1.0, 0.32), fontsize=FS_LEGEND,
    borderaxespad=0.0, handletextpad=0.4)
add_panel_label(axA, "A", x=-0.245, y=1.02)

# ------------------------------ Panel B -----------------------------------
axB = fig.add_subplot(gs[1, 0])
n = len(P)
for r in P:
    col = (C_SIG if r["survives"] else
           (C_DARK if r["q"] < 0.05 else "#C4C4C4"))
    axB.plot([0, r["diff"]], [r["rank"] - 1, r["rank"] - 1],
             color="#E0E0E0" if r["q"] >= 0.05 else "#B9B9B9",
             lw=0.7, zorder=1)
    axB.scatter(r["diff"], r["rank"] - 1, s=17 if r["q"] < 0.05 else 9,
                marker="o", facecolor=col,
                edgecolor=C_DARK if r["survives"] else col, linewidths=0.5,
                zorder=3)
axB.axvline(0, color=C_DARK, lw=0.9, zorder=2)
axB.axhline(len(sig) - 0.5, color=C_DARK, lw=0.8, ls="--", zorder=2)
axB.text(-3.3, 7.4, "BH-FDR < 0.05 (6 tests above the line)", fontsize=6.8,
         color=C_DARK, ha="left")
for ct, ea, eb, ty in [("IPhC_IBC", "CON", "SEN", 4.3),
                       ("HC", "CON", "SEN", 15.0)]:
    r = next(x for x in P if x["cell_type"] == ct and x["ep_A"] == ea
             and x["ep_B"] == eb)
    axB.annotate("%s (%s vs %s)" % (ct, ea, eb),
                 xy=(r["diff"], r["rank"] - 1),
                 xytext=(-3.3, ty), fontsize=6.8,
                 fontweight="bold" if r["q"] < 0.05 else "normal",
                 color=C_DARK, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", lw=0.5, color="#888888",
                                 shrinkA=0, shrinkB=2))
axB.set_ylim(n - 0.5, -0.5)
axB.set_xlim(*LIM)
axB.set_xlabel("D = q95(A) $-$ q95(B)", fontsize=FS_LABEL)
axB.set_ylabel("Pairwise tests ranked by Monte-Carlo p\n(1 = most significant)",
               fontsize=FS_LABEL)
axB.set_yticks([1, 10, 20, 30, 40, 50, 60, 70, 84])
axB.grid(axis="x", ls=":", lw=0.4, color="#EDEDED", zorder=0)
axB.tick_params(axis="both", labelsize=FS_TICK)
add_panel_label(axB, "B", x=-0.185, y=1.02)

png, pdf = save(fig, "Fig2_cross_endpoint")
print("wrote", png)
print("wrote", pdf)
print("n tests =", len(P), "| BH<0.05 =", len(sig))
for r in top:
    print("  %-12s %-5s vs %-5s  D=%+.4f  mc_p=%.4g  q=%.4g  survives=%s"
          % (r["cell_type"], r["ep_A"], r["ep_B"], r["diff"], r["mc_p"],
             r["q"], r["survives"]))
