# -*- coding: utf-8 -*-
"""Shared plotting style for the manuscript figures.

Light theme (white background, dark text), generic sans-serif fonts only
(DejaVu Sans / Arial). No CJK glyphs anywhere inside the figures.
Colour palette follows the colour-blind-safe Okabe-Ito scheme and every
categorical distinction is additionally encoded by marker shape, so colour
is never the sole channel.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- fonts (generic sans-serif only) ------------------------------------
FONT = {"family": "sans-serif",
        "sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "Liberation Sans"],
        "size": 8.0}
matplotlib.rc("font", **FONT)
matplotlib.rc("mathtext", fontset="dejavusans", default="regular")
# NOTE (MI-22): text that must survive PDF text extraction / copy-paste is
# written as plain Unicode (e.g. U+2032 PRIME) rather than mathtext, because
# mathtext glyphs are not mapped back to Unicode in the PDF text layer.
matplotlib.rc("pdf", fonttype=42)          # embed TrueType in PDF
matplotlib.rc("ps", fonttype=42)
matplotlib.rc("svg", fonttype="none")
matplotlib.rc("axes", linewidth=0.8, edgecolor="#333333",
              labelsize=8.0, titlesize=9.0, labelcolor="#1A1A1A",
              titlecolor="#1A1A1A", axisbelow=True)
matplotlib.rc("xtick", labelsize=7.0, color="#1A1A1A", direction="out")
matplotlib.rc("ytick", labelsize=7.0, color="#1A1A1A", direction="out")
matplotlib.rc("xtick.major", width=0.8, size=3.0)
matplotlib.rc("ytick.major", width=0.8, size=3.0)
matplotlib.rc("legend", fontsize=7.0, frameon=False, handlelength=1.2,
              handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
matplotlib.rc("figure", facecolor="white", dpi=400)
matplotlib.rc("savefig", facecolor="white", dpi=400)
matplotlib.rc("lines", linewidth=1.0, markersize=4.0)
matplotlib.rc("grid", color="#E3E3E3", linewidth=0.5)

# ---- palette (Okabe-Ito, colour-blind safe) -----------------------------
C_SEN = "#0072B2"      # blue      - primary endpoint SEN
C_CON = "#E69F00"      # orange    - exploratory endpoint CON
C_NOISE = "#CC79A7"    # purple    - exploratory endpoint NOISE
C_SIG = "#D55E00"      # vermillion - highlight
C_GREY = "#8C8C8C"
C_GREY_LIGHT = "#C9C9C9"
C_DARK = "#1A1A1A"
C_LINE = "#333333"
C_MANH_A = "#A9C0DE"
C_MANH_B = "#5C7FA8"
C_REF = "#B0B0B0"

EP_COLOR = {"SEN": C_SEN, "CON": C_CON, "NOISE": C_NOISE}
EP_MARKER = {"SEN": "o", "CON": "s", "NOISE": "^"}
EP_LABEL = {"SEN": "SEN (sensorineural)", "CON": "CON (conductive, exploratory)",
            "NOISE": "NOISE (noise-induced, exploratory)"}

OUT_DIR = "C:/Users/docta/WorkBuddy/纯生信/figures"

DPI = 400


def save(fig, name):
    """Save one figure as PNG (400 dpi) and PDF (vector)."""
    import os
    os.makedirs(OUT_DIR, exist_ok=True)
    png = os.path.join(OUT_DIR, name + ".png")
    pdf = os.path.join(OUT_DIR, name + ".pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(pdf, dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return png, pdf


def add_panel_label(ax, label, x=-0.02, y=1.04):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="right", color=C_DARK)


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        if s not in keep:
            ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", which="both", length=2.5)
