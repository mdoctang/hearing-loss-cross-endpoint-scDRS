# -*- coding: utf-8 -*-
"""Figure 4 (final numbering; file fig4_sensitivity.py).

Fig 2 - sensitivity of the IPhC_IBC CON-vs-SEN difference.

All panels show D = q95(CON) - q95(SEN) for the inner phalangeal cell /
inner border cell (IPhC_IBC) group in the whole-cochlea atlas.

Layout note (MA-7): the three panels are stacked vertically and the figure is
7.1 in wide (BMC double-column text width), so no scaling is needed at
production size; all text is >= 7 pt and remains >= 7 pt after insertion.

The number of pairwise tests entering the BH correction is NOT the same in
every re-analysis, so it is printed next to every q value (m = number of tests
in the correction family). The top-2000 gene-set analysis could only be
completed for SEN and CON (NOISE endpoint not available at that gene-set
size), so its BH family contains 28 tests instead of 84; that point is drawn
with a hollow marker plus an explicit "not comparable" note, and a
conservative upper bound on the 84-test-family q is printed for it.

Panel A: cell-level re-aggregation under 5 configurations.
Panel B: gene-set size gradient (top 500 / 1000 / 2000 MAGMA genes).
Panel C: gene-set construction variants.

Inputs (read-only):
  project/analysis/scdrs_main/whole_cochlea_atlas/compare_pairwise_tests.csv
  project/analysis/sensitivity/reaggregate_whole_cochlea_atlas/pairwise_*.csv
  project/analysis/sensitivity/gs500/whole_cochlea_atlas/compare_pairwise_tests.csv
  project/analysis/sensitivity/gs2000/whole_cochlea_atlas/PARTIAL_SENCON_only_pairwise_tests.csv
  project/analysis/sensitivity/mhc_excluded/whole_cochlea_atlas/compare_pairwise_tests.csv
  project/analysis/sensitivity/ortholog_builtin/whole_cochlea_atlas/compare_pairwise_tests.csv
"""
import csv, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import C_GREY, C_DARK, save, add_panel_label

ANA = "C:/Users/docta/WorkBuddy/纯生信/project/analysis"
SENS = os.path.join(ANA, "sensitivity")
MAIN = os.path.join(ANA, "scdrs_main", "whole_cochlea_atlas",
                    "compare_pairwise_tests.csv")
FULL_FAMILY = 84


def iphc(path):
    """IPhC_IBC CON-vs-SEN row plus the size of the BH correction family."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        if (r["cell_type"] == "IPhC_IBC" and r["ep_A"] == "CON"
                and r["ep_B"] == "SEN"):
            return dict(D=float(r["diff"]), mc_p=float(r["mc_p"]),
                        q=float(r["qval_BH"]), m=len(rows))
    raise KeyError(path)


# ---- Panel A: cell-level re-aggregation configs ---------------------------
R = os.path.join(SENS, "reaggregate_whole_cochlea_atlas")
A = [("All cells (primary)", os.path.join(R, "pairwise_all.csv")),
     ("Drop Miao2024", os.path.join(R, "pairwise_dataset-not-Miao2024.csv")),
     ("Drop Sun2023_full",
      os.path.join(R, "pairwise_dataset-not-Sun2023_full.csv")),
     ("Control condition only",
      os.path.join(R, "pairwise_condition-Control.csv")),
     ("200 control gene sets", os.path.join(R, "pairwise_ctrl_sub=200.csv"))]
rows_a = [dict(iphc(p), label=lab) for lab, p in A]

# ---- Panel B: gene-set size ----------------------------------------------
B = [("Top 500 genes", os.path.join(SENS, "gs500", "whole_cochlea_atlas",
                                    "compare_pairwise_tests.csv")),
     ("Top 1000 genes\n(primary)", MAIN),
     ("Top 2000 genes", os.path.join(SENS, "gs2000", "whole_cochlea_atlas",
                                     "PARTIAL_SENCON_only_pairwise_tests.csv"))]
rows_b = [dict(iphc(p), label=lab) for lab, p in B]

# ---- Panel C: gene-set construction variants -----------------------------
C = [("Primary gene set", MAIN),
     ("MHC region excluded\n(271 genes dropped)",
      os.path.join(SENS, "mhc_excluded", "whole_cochlea_atlas",
                   "compare_pairwise_tests.csv")),
     ("scDRS built-in\northolog bridge",
      os.path.join(SENS, "ortholog_builtin", "whole_cochlea_atlas",
                   "compare_pairwise_tests.csv"))]
rows_c = [dict(iphc(p), label=lab) for lab, p in C]

# conservative upper bound on the 84-family q for the 28-family point
# (BH q = p x m / rank; the largest value is obtained at the most adverse rank 1)
q84_ub = min(1.0, rows_b[2]["mc_p"] * FULL_FAMILY / 1.0)

# =============================== FIGURE ===================================
# MA-7: single-column-stacked layout, 7.1 in wide = BMC double-column text
# width; every font size >= 7 pt so nothing falls below ~6 pt when printed.
FS_ANNOT = 7.5     # q annotations
FS_NOTE = 7.0      # top-2000 caveat block
FS_TICK = 7.0
FSYLABEL = 8.0
FS_TITLE = 9.0
FS_LEGEND = 7.0

fig = plt.figure(figsize=(6.69, 8.75))
gs = gridspec.GridSpec(3, 1, height_ratios=[1.22, 0.95, 0.95], hspace=0.60,
                       left=0.19, right=0.975, top=0.95, bottom=0.135)
LIM = (-2.7, 1.7)


def draw_panel(ax, rows, title, note_row=None):
    ys = np.arange(len(rows))
    for i, r in enumerate(rows):
        hollow = (note_row is not None and i == note_row)
        col = C_DARK if r["q"] < 0.05 else C_GREY
        ax.plot([LIM[0], r["D"]], [i, i], color="#DDDDDD", lw=1.0, zorder=1)
        ax.scatter(r["D"], i, s=58 if hollow else 48,
                   marker="D" if hollow else "o",
                   facecolor="white" if hollow else col,
                   edgecolor=C_DARK if hollow else col,
                   linewidths=1.3 if hollow else 0.7, zorder=3)
        if hollow:
            txt = ("q = %.4g  (%d-comparison family, SEN/CON only)\n"
                   "not comparable to the %d-comparison q;\n"
                   "conservative upper bound on the %d-family q \u2264 %.1e"
                   % (r["q"], r["m"], FULL_FAMILY, FULL_FAMILY, q84_ub))
            ax.annotate(txt, xy=(-1.70, i), xycoords="data", fontsize=FS_NOTE,
                        color=C_DARK, va="center", ha="left", linespacing=1.45)
        else:
            ax.text(r["D"] + 0.13, i, "q = %.4g  (m = %d tests)"
                    % (r["q"], r["m"]), fontsize=FS_ANNOT,
                    color=C_DARK if r["q"] < 0.05 else "#666666",
                    va="center", ha="left",
                    fontweight="bold" if r["q"] * 2 < 0.05 else "normal")
    ax.axvline(0, color=C_DARK, lw=1.0, zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=FS_TICK)
    ax.set_ylim(len(rows) - 0.45, -0.75)
    ax.set_xlim(*LIM)
    ax.set_title(title, fontsize=FS_TITLE, pad=7)
    ax.grid(axis="x", ls=":", lw=0.5, color="#EDEDED", zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=FS_TICK)


axA = fig.add_subplot(gs[0, 0])
draw_panel(axA, rows_a, "Cell-level re-aggregation")
axA.set_xlabel("D = q95(CON) $-$ q95(SEN)", fontsize=FSYLABEL)
add_panel_label(axA, "A", x=-0.155, y=1.04)

axB = fig.add_subplot(gs[1, 0])
draw_panel(axB, rows_b, "Gene-set size", note_row=2)
add_panel_label(axB, "B", x=-0.155, y=1.06)

axC = fig.add_subplot(gs[2, 0])
draw_panel(axC, rows_c, "Gene-set construction")
axC.set_xlabel("D = q95(CON) $-$ q95(SEN)", fontsize=FSYLABEL)
add_panel_label(axC, "C", x=-0.155, y=1.06)

fig.legend(handles=[
    Line2D([], [], marker="o", ls="none", mfc=C_DARK, mec=C_DARK, ms=6.5,
           label="BH-FDR < 0.05 (within analysis)"),
    Line2D([], [], marker="o", ls="none", mfc=C_GREY, mec=C_GREY, ms=6.5,
           label="BH-FDR \u2265 0.05 (within analysis)"),
    Line2D([], [], marker="D", ls="none", mfc="white", mec=C_DARK, ms=7,
           label="28-comparison family only (SEN/CON); q NOT comparable to the "
                 "84-comparison q")],
    loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=1, fontsize=FS_LEGEND,
    frameon=False, handletextpad=0.4, labelspacing=0.35, borderaxespad=0.0)

png, pdf = save(fig, "Fig4_sensitivity")
print("wrote", png)
print("wrote", pdf)
print("conservative upper bound on 84-family q for gs2000: %.4e" % q84_ub)
for nm, rows in [("A", rows_a), ("B", rows_b), ("C", rows_c)]:
    for r in rows:
        print("  %s  %-34s D=%+.4f  mc_p=%.4g  q=%.4g  m=%d"
              % (nm, r["label"].replace("\n", " "), r["D"], r["mc_p"],
                 r["q"], r["m"]))
