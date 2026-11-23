"""Figure of the full-domain GOAT grid of TCV-X21 and its source regions.

usage (from us_runs/): python tcvx21/full_grid_fig.py <full grid with sources> <out.png>

Left: the whole grid; cells that the standard grid also has in grey, the cells
added inside its core boundary (psiN 0.90) in blue. Middle: where the heating
and the fuelling are deposited. Right: the cells round the magnetic axis.
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, netCDF4, flatgrid
import compare as c
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

BLUE, GREY, INK, MUTED = "#2a78d6", "#b9b8b0", "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e", "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.titlesize": 10.5, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
WALL = np.array([(0.671, 0.750), (0.624, 0.700), (0.624, -0.700), (0.671, -0.750), (0.965, -0.750), (1.136, -0.543), (1.136, 0.543), (0.965, 0.750), (0.671, 0.750)])

path, out = sys.argv[1:3]
g = flatgrid.read(path); n = g["nreal"]
polys = [np.stack([g["vertex_R"][v[v >= 0]], g["vertex_Z"][v[v >= 0]]], 1) for v in g["cell_vertices"]][:n]
psin = 1 + g["psi"][:n]*2*np.pi/c.q["simagx"]
tube = g["tube"][:n]
closed = np.array([(g["prev"][:n][tube == t] >= 0).all() for t in range(tube.max() + 1)])[tube]
new = closed & (psin < 0.9)
with netCDF4.Dataset(path) as d:
    s_of_l = np.asarray(d["us_logical"][:])
    p_src = np.asarray(d["Pd+_src"][0, :])[s_of_l[:n]]
    n_src = np.asarray(d["Nd+_src"][0, :])[s_of_l[:n]]
RR, ZZ = np.meshgrid(np.linspace(c.er[0], c.er[-1], 300), np.linspace(c.ez[0], c.ez[-1], 600), indexing="ij")
PSN = c.psin_RZ(RR, ZZ)

fig, ax = plt.subplots(1, 3, figsize=(13, 7.2), gridspec_kw={"width_ratios": [1, 1, 1.25]})
for a in ax:
    a.plot(WALL[:, 0], WALL[:, 1], color=INK, lw=1.2, zorder=5)
    a.set_aspect("equal"); a.set_xlabel("R [m]")
a = ax[0]
a.add_collection(PolyCollection([p for p, m in zip(polys, new) if not m], facecolors="none", edgecolors=GREY, linewidths=0.3))
a.add_collection(PolyCollection([p for p, m in zip(polys, new) if m], facecolors="none", edgecolors=BLUE, linewidths=0.45))
a.contour(RR, ZZ, np.where(ZZ > -0.42, PSN, np.nan), levels=[0.9], colors=[INK], linewidths=0.9, linestyles="--", zorder=4)
a.set_xlim(0.59, 1.16); a.set_ylim(-0.78, 0.78); a.set_ylabel("Z [m]")
a.set_title("Whole-plasma grid, %d cells" % n)
a = ax[1]
heat = p_src > 0; fuel = n_src > 0
a.add_collection(PolyCollection([p for p, m in zip(polys, heat) if m], facecolors=BLUE, alpha=0.55, edgecolors="none"))
a.add_collection(PolyCollection([p for p, m in zip(polys, fuel) if m], facecolors="#eb6834", alpha=0.7, edgecolors="none"))
a.add_collection(PolyCollection(polys, facecolors="none", edgecolors=GREY, linewidths=0.2))
a.contour(RR, ZZ, PSN, levels=[1.0], colors=[INK], linewidths=0.7, linestyles="--", zorder=4)
a.set_xlim(0.59, 1.16); a.set_ylim(-0.78, 0.78); a.set_yticklabels([])
a.set_title("Source regions")
a = ax[2]
a.add_collection(PolyCollection(polys, facecolors="none", edgecolors=BLUE, linewidths=0.6))
a.plot([c.q["rmagx"]], [c.q["zmagx"]], "+", color=INK, ms=9, mew=1.4)
a.set_xlim(c.q["rmagx"] - 0.08, c.q["rmagx"] + 0.08); a.set_ylim(c.q["zmagx"] - 0.12, c.q["zmagx"] + 0.12)
a.set_title("Cells round the magnetic axis (+)")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)
