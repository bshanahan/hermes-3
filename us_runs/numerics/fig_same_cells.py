"""Figure: the structured code and the flat-cell code on identical cells.

usage (from us_runs/): python numerics/fig_same_cells.py <structured run> <flat run> <flat grid> <out stem>

The orthogonal TCV-X21 grid (32 x 104 cells) is run by the structured
operators, and recast as a list of flat cells (with the structured metric)
for the flat-cell operators. Same inputs: no drifts, no skew correction,
fixed core density and temperature, recycling 0.99, no pumping.
(a) cumulative distribution of the relative difference per cell
(b) outer midplane, (c) outer target: electron density and temperature
"""
import sys
from pathlib import Path

import netCDF4
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, "../tools/unstructured")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as ps  # noqa: E402
import usdata  # noqa: E402

ps.setup()
plt = ps.plt
spath, fpath, gridfile, stem = sys.argv[1:5]


def structured(name, it=-1):
    out = []
    files = sorted(Path(spath).glob("BOUT.dmp.*.nc"), key=lambda p: int(p.name.split(".")[2]))
    # assemble the global (x, y) array from the processors' pieces
    with netCDF4.Dataset(str(files[0])) as d:
        nx, ny = int(d["nx"][...]), int(d["ny"][...])
        mxg, myg = int(d["MXG"][...]), int(d["MYG"][...])
        nxpe = int(d["NXPE"][...])
    a = np.full((nx, ny), np.nan)
    for p in files:
        with netCDF4.Dataset(str(p)) as d:
            v = d[name]
            loc = np.asarray(v[it, :, :, 0])
            if "conversion" in v.ncattrs():
                loc = loc * v.getncattr("conversion")
            mxsub, mysub = int(d["MXSUB"][...]), int(d["MYSUB"][...])
            px, py = int(d["PE_XIND"][...]), int(d["PE_YIND"][...])
            x0 = px * mxsub
            y0 = py * mysub
            xs = 0 if px == 0 else mxg
            xe = mxsub + 2 * mxg if px == nxpe - 1 else mxsub + mxg
            gx0 = x0 + xs
            a[gx0:gx0 + (xe - xs), y0:y0 + mysub] = loc[xs:xe, myg:myg + mysub]
    return a


r = usdata.Run(fpath, gridfile)
n = r.n
import flatgrid  # noqa: E402
_g = flatgrid.read(gridfile)
xi = np.asarray(_g["structured_x"][:n]).astype(int)      # x index including the guard cells
yj = np.asarray(_g["structured_y"][:n]).astype(int)
qe = 1.602176634e-19

fields = [("Nd+", r"ion density", ps.BLUE, "-"), ("Te", r"electron temperature", ps.ORANGE, "--"),
          ("Td+", r"ion temperature", ps.AQUA, "-."), ("Nd", r"atom density", ps.YELLOW, ":")]
vals = {}
for name, label, color, ls in fields:
    s = structured(name)[xi, yj]
    u = np.asarray(r(name))[-1][:n]
    vals[name] = (s, u)

fig = plt.figure(figsize=(ps.WIDTH, 3.9))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], hspace=0.55, wspace=0.42)
ax = fig.add_subplot(gs[:, 0])
for name, label, color, ls in fields:
    s, u = vals[name]
    diff = np.sort(np.abs(u / s - 1))
    frac = np.arange(1, len(diff) + 1) / len(diff)
    ax.plot(np.maximum(diff, 1e-7), frac, color=color, ls=ls, lw=1.4, label=label)
    print("%-22s median %.2e  90%% %.2e  max %.2e" % (label, np.median(diff), np.percentile(diff, 90), diff.max()))
ax.set_xscale("log")
ax.set_xlim(1e-5, 1)
ax.set_ylim(0, 1.02)
ax.set_xlabel(r"relative difference per cell")
ax.set_ylabel(r"fraction of cells")
ax.legend(loc="upper left", fontsize=7.5, handlelength=2.0)
ax.set_title(r"(a) All %d cells" % n, loc="left", fontsize=9, pad=4)

# outer midplane and outer target: the structured row through the outer
# midplane (largest R near Z = 0) and the last row of the outer leg
g = r.g
R, Z = g["R"][:n], g["Z"][:n]
ny = yj.max() + 1
cand = np.flatnonzero(np.abs(Z) < 0.03)
j_omp = yj[cand[np.argmax(R[cand])]]
panels = [(0, 1, j_omp, "Nd+", 1e-19, r"$n_e$ [$10^{19}\,$m$^{-3}$]", r"(b) Outer midplane"),
          (1, 1, j_omp, "Te", 1.0, r"$T_e$ [eV]", r"(c) Outer midplane"),
          (0, 2, ny - 1, "Nd+", 1e-19, r"$n_e$ [$10^{19}\,$m$^{-3}$]", r"(d) Outer target"),
          (1, 2, ny - 1, "Te", 1.0, r"$T_e$ [eV]", r"(e) Outer target")]
for row, col, j, name, scale, ylabel, title in panels:
    ax = fig.add_subplot(gs[row, col])
    sel = np.flatnonzero(yj == j)
    sel = sel[np.argsort(xi[sel])]
    x = R[sel] if "midplane" in title else Z[sel]
    color = ps.BLUE if name == "Nd+" else ps.ORANGE
    ax.plot(x, vals[name][0][sel] * scale, color=color, lw=1.4)
    ax.plot(x, vals[name][1][sel] * scale, ls="none", marker="o", ms=3.0, mfc="white", mec=color, mew=0.9)
    ax.set_ylabel(ylabel)
    if row == 1:
        ax.set_xlabel(r"$R$ [m]" if "midplane" in title else r"$Z$ [m]")
    ax.set_title(title, loc="left", fontsize=9, pad=4)
    ax.set_ylim(bottom=0)
handles = [Line2D([], [], color=ps.SECONDARY, lw=1.4, label=r"structured code"),
           Line2D([], [], color=ps.SECONDARY, ls="none", marker="o", ms=3.0, mfc="white", mec=ps.SECONDARY,
                  label=r"flat-cell code, same cells")]
fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.62, 1.04))
ps.save(fig, stem)
print("wrote", stem)
