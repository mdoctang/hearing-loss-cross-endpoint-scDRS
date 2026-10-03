# -*- coding: utf-8 -*-
"""Figure 1 (final numbering; file fig1_gene_level.py).

Fig 1 - gene-level MAGMA results for the three FinnGen R13 hearing-loss endpoints.

Layout note (BMC compliance): single-column figure, plotted at 170 mm wide
(BMC full-page width) so that nominal font size == effective font size when
placed in the journal; all text >= 7 pt; height kept well inside the 225 mm
figure-and-legend cap.

Panel A-C: three stacked Manhattan plots of gene-level -log10(P) (MAGMA gene
          test, FUMA output, GRCh37 coordinates), shared x and y axes.
Panel D:   gene-level z-statistic distributions with the N(0,1) reference.

Inputs (read-only):
  project/analysis/fuma_output/<endpoint>/<endpoint>.magma.genes.out
"""
import csv, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (C_SIG, C_DARK, C_MANH_A, C_MANH_B, save, add_panel_label,
                      check_text_collisions)

# every text artist that must not overlap another (checked before saving)
TEXT_ITEMS = []

BASE = "C:/Users/docta/WorkBuddy/纯生信/project/analysis/fuma_output"
BONF = 2.634352e-06                      # Bonferroni 0.05 / 18,980 genes
BONF_LOG = -np.log10(BONF)

ENDPOINTS = [("SEN", "H8_HL_SEN_NAS"), ("CON", "H8_HL_CON_NAS"),
             ("NOISE", "H8_NOISEINNER")]

ANNO_SEN = ["ARHGEF28", "NOL12", "HLA-DQA1", "EML6"]
ANNO_CON = ["ACOXL", "PARP12"]
ANNO_NOISE = ["UBE3A"]


def load(endpoint_dir):
    path = os.path.join(BASE, endpoint_dir, endpoint_dir + ".magma.genes.out")
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return dict(
        chrom=np.array([int(r["CHR"]) for r in rows]),
        start=np.array([int(r["START"]) for r in rows]),
        stop=np.array([int(r["STOP"]) for r in rows]),
        p=np.array([float(r["P"]) for r in rows]),
        z=np.array([float(r["ZSTAT"]) for r in rows]),
        sym=np.array([r["SYMBOL"] for r in rows]))


data = {ep: load(d) for ep, d in ENDPOINTS}

chroms = list(range(1, 23))
d0 = data["SEN"]
chr_len = {c: int(d0["stop"][d0["chrom"] == c].max()) for c in chroms}
offset, cum = {}, 0
for c in chroms:
    offset[c] = cum
    cum += chr_len[c]
GENOME_LEN = cum
chr_mid = {c: offset[c] + chr_len[c] / 2.0 for c in chroms}


def to_genomic(chrom, start):
    return np.array([offset[c] + s for c, s in zip(chrom, start)], dtype=float)


# =============================== FIGURE ===================================
FS_TICK = 7.0
FS_SMALL = 7.0
FS_LABEL = 8.0
FS_BOX = 8.0

fig = plt.figure(figsize=(6.69, 7.3))
gs = gridspec.GridSpec(4, 1, height_ratios=[1.0, 1.0, 1.0, 1.35], hspace=0.46,
                       left=0.088, right=0.985, top=0.935, bottom=0.072)

panel_axes = []
for i, (ep, _) in enumerate(ENDPOINTS):
    ax = fig.add_subplot(gs[i, 0],
                         sharex=(panel_axes[0] if panel_axes else None),
                         sharey=(panel_axes[0] if panel_axes else None))
    panel_axes.append(ax)
    d = data[ep]
    x = to_genomic(d["chrom"], d["start"])
    y = -np.log10(d["p"])
    sig = d["p"] < BONF
    n_sig = int(sig.sum())

    for k, c in enumerate(chroms):
        m = d["chrom"] == c
        ax.scatter(x[m][~sig[m]], y[m][~sig[m]], s=2.2,
                   c=(C_MANH_A if k % 2 == 0 else C_MANH_B),
                   linewidths=0, marker="o", zorder=2)
    if sig.any():
        ax.scatter(x[sig], y[sig], s=11.0, c=C_SIG, linewidths=0, marker="o",
                   zorder=3)

    ax.axhline(BONF_LOG, color=C_DARK, lw=0.9, ls="--", zorder=4)
    ax.set_xlim(0, GENOME_LEN)
    ax.set_ylim(-0.4, 23.4)
    ax.set_ylabel("-\u2009log$_{10}$(P)\ngene-level MAGMA", fontsize=FS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.grid(axis="y", ls=":", lw=0.4, color="#EDEDED", zorder=0)

    top_p = d["p"].min()
    top_sym = d["sym"][int(np.argmin(d["p"]))]
    tag = ("%s — %d / %s genes pass Bonferroni" % (ep, n_sig, f"{len(d['p']):,}")
           if n_sig else "%s — no gene passes Bonferroni" % ep)
    tag_txt = ax.text(0.012, 0.90, tag, transform=ax.transAxes, fontsize=FS_BOX,
                      fontweight="bold", color=C_DARK, va="top",
                      bbox=dict(boxstyle="round,pad=0.3", fc="white",
                                ec="#CCCCCC", lw=0.6))
    TEXT_ITEMS.append((ep + " tag box", tag_txt))
    # The min-P line lives in the top-RIGHT corner, not under the tag box.
    # Reason (2026-10-03): at 0.68 axes height, left-aligned, it overlapped the
    # EML6 gene label; the top-right corner of every panel is empty because the
    # only far-right annotation (NOL12, panel A) is placed below its point.
    minp_txt = ax.text(0.985, 0.90, "min gene P = %.3g (%s)" % (top_p, top_sym),
                       transform=ax.transAxes, fontsize=FS_SMALL,
                       color="#555555", va="top", ha="right")
    TEXT_ITEMS.append((ep + " minP line", minp_txt))

    # Genes whose label must go BELOW the point.  ARHGEF28 stays above (it is
    # the single topmost gene); the rest sit under the tag box, so a label
    # placed above them would land inside the box rectangle
    # (axes x 0.012-0.484, y 0.718-0.900) and collide with the bold summary.
    BELOW = {"EML6", "HLA-DQA1", "NOL12"}
    for j, name in enumerate({"SEN": ANNO_SEN, "CON": ANNO_CON,
                              "NOISE": ANNO_NOISE}[ep]):
        idx = np.where(d["sym"] == name)[0]
        if idx.size == 0:
            continue
        idx = idx[0]
        gx, gy = float(x[idx]), float(y[idx])
        ha = "center"
        if gx > 0.93 * GENOME_LEN:
            ha = "right"
            gx = GENOME_LEN * 0.997
        if name in BELOW:
            ann = ax.annotate(name, xy=(gx, gy), xytext=(gx, gy - 1.2),
                              fontsize=FS_SMALL, color=C_DARK, ha=ha, va="top",
                              arrowprops=dict(arrowstyle="-", lw=0.5,
                                              color="#888888", shrinkA=0,
                                              shrinkB=1.5))
        else:
            ann = ax.annotate(name, xy=(gx, gy),
                              xytext=(gx, gy + 1.0 + 0.85 * (j % 3)),
                              fontsize=FS_SMALL, color=C_DARK, ha=ha,
                              va="bottom",
                              arrowprops=dict(arrowstyle="-", lw=0.5,
                                              color="#888888", shrinkA=0,
                                              shrinkB=1.5))
        TEXT_ITEMS.append((ep + " " + name, ann))
    add_panel_label(ax, ["a", "b", "c"][i], x=-0.008, y=1.015)

axm = panel_axes[-1]
axm.set_xticks([chr_mid[c] for c in chroms if c % 2 == 1])
axm.set_xticklabels([str(c) for c in chroms if c % 2 == 1], fontsize=FS_TICK)
axm.set_xlabel("Chromosome (GRCh37 gene coordinates)", fontsize=FS_LABEL)
axm.tick_params(axis="x", length=0)
thr_txt = panel_axes[1].text(0.988, 0.62, "Bonferroni threshold P = 2.63e-6",
                             transform=panel_axes[1].transAxes,
                             fontsize=FS_SMALL, color=C_DARK, ha="right",
                             va="bottom")
TEXT_ITEMS.append(("threshold note", thr_txt))

# ---- Panel D: gene-level z distribution ---------------------------------
axd = fig.add_subplot(gs[3, 0])
bins = np.linspace(-6, 8, 141)
for ep, colour in [("SEN", "#0072B2"), ("CON", "#E69F00"),
                   ("NOISE", "#CC79A7")]:
    z = data[ep]["z"]
    axd.hist(z, bins=bins, histtype="step", color=colour, lw=1.3,
             label="%s: SD = %.2f, max = %.2f"
                   % (ep, z.std(ddof=1), z.max()),
             density=True, zorder=3)
xs = np.linspace(-6, 8, 600)
axd.plot(xs, np.exp(-xs ** 2 / 2) / np.sqrt(2 * np.pi), color="#B0B0B0",
         lw=1.2, ls="--", zorder=2, label="N(0, 1) reference")
axd.set_xlabel("MAGMA gene-level z statistic", fontsize=FS_LABEL)
axd.set_ylabel("Density", fontsize=FS_LABEL)
axd.set_xlim(-6, 8)
axd.legend(loc="upper left", fontsize=FS_SMALL, frameon=False)
n_txt = axd.text(0.99, 0.97, "n = 18,980 genes per endpoint",
                 transform=axd.transAxes, fontsize=FS_SMALL, color="#555555",
                 ha="right", va="top")
TEXT_ITEMS.append(("panel D n note", n_txt))
axd.grid(axis="y", ls=":", lw=0.4, color="#EDEDED", zorder=0)
axd.tick_params(axis="both", labelsize=FS_TICK)
add_panel_label(axd, "d", x=-0.008, y=1.015)

# ---- hard gate: no two texts may overlap in the rendered figure ----------
check_text_collisions(fig, TEXT_ITEMS)

png, pdf = save(fig, "Fig1_gene_level")
print("wrote", png)
print("wrote", pdf)
for ep, _ in ENDPOINTS:
    d = data[ep]
    print("  %-6s n=%d  n_sig=%d  minP=%.4g (%s)"
          % (ep, len(d["p"]), int((d["p"] < BONF).sum()), d["p"].min(),
             d["sym"][int(np.argmin(d["p"]))]))
