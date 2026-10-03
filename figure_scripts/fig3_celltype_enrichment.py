# -*- coding: utf-8 -*-
"""Figure 3 (final numbering; file fig2_celltype_enrichment.py).

Fig 4 - cell-type-level scDRS enrichment in the whole-cochlea atlas.

Layout note (BMC compliance): single-column figure, plotted at 170 mm wide so
that nominal font size == effective font size; all text >= 7 pt; height inside
the 225 mm cap. The two panels are stacked vertically.

Panel A: q95 of the normalized scDRS score for 28 cell types x 3 endpoints,
         sorted by q95_SEN. Filled marker = endpoint-level enrichment at
         nominal MC p < 0.05 (within-endpoint test); double dagger = survives
         the prior two-family correction (q' = 2q < 0.05).
Panel B: q95(SEN) against q95(CON) with the identity line. Points above the
         line have D = q95(CON) - q95(SEN) < 0, i.e. SEN > CON.

Input (read-only):
  project/analysis/scdrs_main/whole_cochlea_atlas/compare_endpoint_matrix.csv
  project/analysis/scdrs_main/whole_cochlea_atlas/compare_pairwise_tests.csv
"""
import csv, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (C_SEN, C_CON, C_NOISE, C_SIG, C_DARK, EP_COLOR,
                      EP_MARKER, save, add_panel_label)

ATLAS = ("C:/Users/docta/WorkBuddy/纯生信/project/analysis/scdrs_main/"
         "whole_cochlea_atlas")

with open(os.path.join(ATLAS, "compare_endpoint_matrix.csv"), newline="") as fh:
    M = list(csv.DictReader(fh))
with open(os.path.join(ATLAS, "compare_pairwise_tests.csv"), newline="") as fh:
    P = list(csv.DictReader(fh))

names = [r["cell_type"] for r in M]
n_cell = {r["cell_type"]: int(float(r["n_cell"])) for r in M}
q95 = {ep: {r["cell_type"]: float(r["q95_" + ep]) for r in M}
       for ep in ("SEN", "CON", "NOISE")}
mcp = {ep: {r["cell_type"]: float(r["mcp_" + ep]) for r in M}
       for ep in ("SEN", "CON", "NOISE")}

bh005 = {r["cell_type"] for r in P if float(r["qval_BH"]) < 0.05}
prior = {r["cell_type"] for r in P if float(r["qval_BH"]) * 2 < 0.05}

order = sorted(names, key=lambda c: -q95["SEN"][c])

FS_TICK = 7.0
FS_ANNOT = 7.0
FS_LABEL = 8.0
FS_LEGEND = 7.0

# =============================== FIGURE ===================================
fig = plt.figure(figsize=(6.69, 8.95))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.42, 1.0], hspace=0.88,
                       left=0.185, right=0.985, top=0.95, bottom=0.058)

# ------------------------------ Panel A -----------------------------------
axA = fig.add_subplot(gs[0, 0])
ys = np.arange(len(order))
for i, ct in enumerate(order):
    vals = [q95["CON"][ct], q95["NOISE"][ct], q95["SEN"][ct]]
    axA.plot([min(vals), max(vals)], [i, i], color="#DDDDDD", lw=0.9, zorder=1)
    for ep in ("SEN", "CON", "NOISE"):
        filled = mcp[ep][ct] < 0.05
        axA.scatter(q95[ep][ct], i, s=26 if filled else 21,
                    marker=EP_MARKER[ep],
                    facecolor=EP_COLOR[ep] if filled else "white",
                    edgecolor=EP_COLOR[ep] if filled else "#9A9A9A",
                    linewidths=0.8, zorder=3)

axA.set_yticks(ys)
axA.set_yticklabels(
    [(ct + " (n=%s) \u2021" % f"{n_cell[ct]:,}") if ct in prior else
     (ct + " (n=%s)" % f"{n_cell[ct]:,}") for ct in order],
    fontsize=FS_TICK)
for tick, ct in zip(axA.get_yticklabels(), order):
    tick.set_fontweight("bold" if ct in bh005 else "normal")
axA.invert_yaxis()
axA.set_xlabel("95th percentile of normalized scDRS score (q95)",
               fontsize=FS_LABEL)
axA.set_xlim(0.4, 4.45)
axA.grid(axis="x", ls=":", lw=0.4, color="#EDEDED", zorder=0)
axA.tick_params(axis="y", length=0)
axA.tick_params(axis="x", labelsize=FS_TICK)

for ct in order:
    if ct in prior:
        axA.scatter(4.36, order.index(ct), marker="*", s=72, color=C_SIG,
                    edgecolor=C_DARK, linewidths=0.5, zorder=4, clip_on=False)

handles = [
    Line2D([], [], marker=EP_MARKER["SEN"], ls="none", mfc=C_SEN, mec=C_SEN,
           ms=5.5, label="SEN (primary endpoint)"),
    Line2D([], [], marker=EP_MARKER["CON"], ls="none", mfc=C_CON, mec=C_CON,
           ms=5.5, label="CON (exploratory)"),
    Line2D([], [], marker=EP_MARKER["NOISE"], ls="none", mfc=C_NOISE,
           mec=C_NOISE, ms=5.5, label="NOISE (exploratory)"),
    Line2D([], [], marker="o", ls="none", mfc="white", mec="#9A9A9A", ms=5.5,
           label="nominal within-endpoint MC p \u2265 0.05 (open symbol)"),
    Line2D([], [], marker="*", ls="none", mfc=C_SIG, mec=C_DARK, ms=8,
           label="survives prior correction (q\u2032 = 2q < 0.05)"),
]
axA.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.33), ncol=2,
           fontsize=FS_LEGEND, frameon=False,
           title="bold: cross-endpoint BH-FDR < 0.05\n"
                 "\u2021: survives prior correction",
           title_fontsize=FS_LEGEND, borderaxespad=0.2, labelspacing=0.3)
add_panel_label(axA, "a", x=-0.125, y=1.015)

# ------------------------------ Panel B -----------------------------------
axB = fig.add_subplot(gs[1, 0])
lim = (0.5, 4.6)
axB.fill_between(lim, lim, [lim[1] * 2] * 2, color="#F5F5F5", zorder=0)
axB.plot(lim, lim, color="#AAAAAA", lw=0.9, ls="--", zorder=1)
key = ["HC", "IPhC_IBC", "DC_PC", "SGC", "Neu", "FC4"]
for ct in names:
    in_bh = ct in bh005
    axB.scatter(q95["CON"][ct], q95["SEN"][ct],
                s=42 if in_bh else 26,
                marker="D" if in_bh else "o",
                facecolor=C_DARK if in_bh else "#BBBBBB",
                edgecolor=C_DARK if in_bh else "#888888",
                linewidths=0.7, zorder=3)
offsets = {"HC": (0.06, -0.05, "left"), "IPhC_IBC": (0.08, 0.02, "left"),
           "DC_PC": (0.08, 0.04, "left"), "SGC": (0.08, -0.10, "left"),
           "Neu": (0.08, -0.02, "left"), "FC4": (0.08, 0.06, "left")}
for ct in key:
    dx, dy, ha = offsets[ct]
    axB.annotate(ct, xy=(q95["CON"][ct], q95["SEN"][ct]),
                 xytext=(q95["CON"][ct] + dx, q95["SEN"][ct] + dy),
                 fontsize=FS_ANNOT, fontweight="bold" if ct in bh005
                 else "normal", color=C_DARK, ha=ha, va="center", zorder=4)
axB.set_xlim(lim)
axB.set_ylim(lim)
axB.set_xlabel("q95(CON), exploratory endpoint", fontsize=FS_LABEL)
axB.set_ylabel("q95(SEN), primary endpoint", fontsize=FS_LABEL)
axB.text(0.045, 0.93, "above the diagonal:\nD = q95(CON) $-$ q95(SEN) $<$ 0",
         transform=axB.transAxes, fontsize=6.8, color="#555555", va="top")
axB.text(0.97, 0.06, "below the diagonal: SEN $<$ CON",
         transform=axB.transAxes, fontsize=6.8, color="#555555",
         ha="right", va="bottom")
axB.grid(ls=":", lw=0.4, color="#EDEDED", zorder=0)
axB.legend(handles=[
    Line2D([], [], marker="D", ls="none", mfc=C_DARK, mec=C_DARK, ms=5.5,
           label="BH-FDR < 0.05 (cross-endpoint)"),
    Line2D([], [], marker="o", ls="none", mfc="#BBBBBB", mec="#888888",
           ms=5, label="not significant")],
    loc="upper right", fontsize=FS_LEGEND, borderaxespad=0.2)
axB.tick_params(axis="both", labelsize=FS_TICK)
add_panel_label(axB, "b", x=-0.125, y=1.015)

png, pdf = save(fig, "Fig3_celltype_enrichment")
print("wrote", png)
print("wrote", pdf)
print("n types:", len(order), "| BH<0.05:", sorted(bh005), "| prior:",
      sorted(prior))
print("top5 by q95_SEN:", [(c, round(q95['SEN'][c], 3)) for c in order[:5]])
