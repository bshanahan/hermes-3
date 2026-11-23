"""Shared figure style for the paper: TeX typesetting in Computer Modern (the
font MathJax uses for TeX), thin marks, categorical colours in a fixed order
(validated for colour-blind separation), each series also told apart by
marker and line style."""
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, SECONDARY, MUTED, RULE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
WIDTH = 7.0      # inches, full text width


def setup():
    mpl.rcParams.update({
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.edgecolor": SECONDARY,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.6,
        "xtick.color": SECONDARY,
        "ytick.color": SECONDARY,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "lines.linewidth": 1.3,
        "lines.markersize": 4.5,
        "lines.markeredgewidth": 0.0,
        "legend.frameon": False,
        "legend.handlelength": 2.4,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def panel_label(ax, text, x=0.03, y=0.95):
    ax.text(x, y, text, transform=ax.transAxes, ha="left", va="top", fontsize=9, color=INK)


def slope_triangle(ax, x0, y0, order, width=0.45, label=None, flip=False):
    """Reference slope on log axes: a short line of the given order from (x0, y0),
    spanning `width` decades in x, with its order written beside it."""
    import numpy as np
    x1 = x0 * 10**width
    y1 = y0 * (x1 / x0)**order
    ax.plot([x0, x1], [y0, y1], color=MUTED, lw=0.8, ls="-", zorder=1)
    xm = np.sqrt(x0 * x1)
    ym = np.sqrt(y0 * y1)
    ax.text(xm, ym * (0.62 if not flip else 1.55), label or r"$\propto h^{%d}$" % order,
            fontsize=7.5, color=SECONDARY, ha="center", va="center")


def save(fig, stem):
    """Save as PDF (for the paper) and PNG (for viewing)"""
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".png")
