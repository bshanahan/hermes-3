"""Figures for the TCV-X21 section of the report.

usage (from us_runs/): python tcvx21/report_figs.py <GOAT run dir> <structured run dir> <structured grid file> <out dir>

Writes grids.png, fields.png, upstream.png, outer_target.png, inner_target.png.
"""
import sys, glob
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, netCDF4, usdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import compare as c

goat_run, struct_run, struct_grid, out = sys.argv[1:5]
BLUE, ORANGE = "#2a78d6", "#eb6834"                    # GOAT grid, orthogonal grid
INK, MUTED, GRIDLINE = "#0b0b0b", "#898781", "#e1e0d9"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e", "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.titlesize": 10.5, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
WALL = np.array([(0.671, 0.750), (0.624, 0.700), (0.624, -0.700), (0.671, -0.750), (0.965, -0.750), (1.136, -0.543), (1.136, 0.543), (0.965, 0.750), (0.671, 0.750)])

# ---- GOAT grid and run
d = c.load(goat_run); g, n = d["g"], d["n"]
goat_polys = [np.stack([g["vertex_R"][v[v >= 0]], g["vertex_Z"][v[v >= 0]]], 1) for v in g["cell_vertices"]][:n]
goat_nn = d["r"]("Nd", True)[-1][:n]

# ---- structured grid and run
with netCDF4.Dataset(struct_grid) as gf:
    cor = [np.stack([np.asarray(gf["Rxy" + s][:]), np.asarray(gf["Zxy" + s][:])], -1)
           for s in ("_corners", "_lower_right_corners", "_upper_right_corners", "_upper_left_corners")]
P = np.stack(cor, 2)                                   # (nx, ny, 4, 2)
files = sorted(glob.glob(struct_run + "/BOUT.dmp.*.nc"), key=lambda f: int(f.split(".")[-2]))
def sget(name):
    parts = []
    for f in files:
        with netCDF4.Dataset(f) as dd:
            v = dd[name]; myg = int(dd["MYG"][...])
            parts.append(np.asarray(v[-1, :, myg:-myg, 0])*(v.getncattr("conversion") if "conversion" in v.ncattrs() else 1.0))
    return np.concatenate(parts, 1)
mxg = 2
s_polys = P[mxg:-mxg].reshape(-1, 4, 2)
s_ne, s_te, s_nn = (sget(k)[mxg:-mxg].ravel() for k in ("Nd+", "Te", "Nd"))

# separatrix from the equilibrium
RR, ZZ = np.meshgrid(np.linspace(c.er[0], c.er[-1], 300), np.linspace(c.ez[0], c.ez[-1], 600), indexing="ij")
PSN = c.psin_RZ(RR, ZZ)

def frame(a, zoom=None, sep=True, lw=1.0):
    a.plot(WALL[:, 0], WALL[:, 1], color=INK, lw=1.3, zorder=5)
    if sep:
        a.contour(RR, ZZ, PSN, levels=[1.0], colors=[INK], linewidths=lw, linestyles="--", zorder=4)
    a.set_aspect("equal")
    if zoom:
        a.set_xlim(zoom[0], zoom[1]); a.set_ylim(zoom[2], zoom[3])
    else:
        a.set_xlim(0.59, 1.16); a.set_ylim(-0.78, 0.78)
    a.set_xlabel("R [m]")

# ---- 1. grids
fig = plt.figure(figsize=(12.5, 8.2))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.45], wspace=0.22, hspace=0.3)
ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[:, 1]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])]
ZOOM = (0.59, 1.01, -0.78, -0.31)
for a, polys, col, title, zoom in ((ax[0], goat_polys, BLUE, "GOAT grid\n%d cells, reaches the wall" % n, None),
                                   (ax[1], s_polys, ORANGE, "Orthogonal grid\n%d cells, ends on flux surfaces" % len(s_polys), None),
                                   (ax[2], goat_polys, BLUE, "GOAT grid, divertor", ZOOM), (ax[3], s_polys, ORANGE, "Orthogonal grid, divertor", ZOOM)):
    a.add_collection(PolyCollection(polys, facecolors="none", edgecolors=col, linewidths=0.35 if zoom is None else 0.6))
    frame(a, zoom, sep=False); a.set_title(title)
    a.set_ylabel("Z [m]")
fig.savefig(out + "/grids.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ---- 2. fields: GOAT | orthogonal for three quantities
fig, ax = plt.subplots(1, 6, figsize=(15.5, 6.6))
quant = [("electron density [m$^{-3}$]", d["Ne"][:n], s_ne, LogNorm(1e16, 2e19), "viridis"),
         ("electron temperature [eV]", d["Te"][:n], s_te, LogNorm(1, 100), "inferno"),
         ("atom density [m$^{-3}$]", goat_nn, s_nn, LogNorm(1e15, 1e19), "viridis")]
for k, (lab, vg, vs, norm, cmap) in enumerate(quant):
    for a, polys, val, name in ((ax[2*k], goat_polys, vg, "GOAT"), (ax[2*k + 1], s_polys, vs, "orthogonal")):
        pc = PolyCollection(polys, array=np.clip(val, norm.vmin, None), norm=norm, cmap=cmap, edgecolors="face", linewidths=0.1)
        a.add_collection(pc); frame(a, lw=0.7); a.set_title(name, fontsize=10)
        if a is not ax[0]: a.set_yticklabels([])
    cb = fig.colorbar(pc, ax=[ax[2*k], ax[2*k + 1]], orientation="horizontal", fraction=0.045, pad=0.09, aspect=30); cb.set_label(lab); cb.outline.set_visible(False)
ax[0].set_ylabel("Z [m]")
fig.savefig(out + "/fields.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ---- 3-5. profiles against the measurements
ds = c.load_structured(struct_run, struct_grid)
runs = (("Hermes-3, GOAT grid", d, BLUE), ("Hermes-3, orthogonal grid", ds, ORANGE))
def profiles(panels, path, ncols, title):
    nrows = int(np.ceil(len(panels)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6*ncols, 3.5*nrows), squeeze=False)
    for a in axes.ravel()[len(panels):]: a.set_visible(False)
    for a, (diag, obs, lab, scale) in zip(axes.ravel(), panels):
        xmin, xmax = np.inf, -np.inf
        for field, col, mk in (("forward", "#52514e", "o"), ("reversed", "#a3a29b", "s")):
            xe, ve, ee = c.exp_profile(field, diag, obs)
            a.errorbar(xe, ve*scale, ee*scale, fmt=mk, ms=3, color=col, lw=0.7, capsize=1.5, label="measured, %s field" % field, zorder=2)
            xmin, xmax = min(xmin, xe.min()), max(xmax, xe.max())
        for name, dd, col in runs:
            xs, vs = c.sim_profile(dd, diag, obs)
            a.plot(xs, vs*scale, "-", color=col, lw=2, label=name, zorder=3)
        pad = 0.12*(xmax - xmin); a.set_xlim(xmin - pad, xmax + pad); a.set_ylim(0, None)
        a.axvline(0, color="#c3c2b7", lw=0.8, zorder=1); a.grid(axis="y", color=GRIDLINE, lw=0.6); a.set_axisbelow(True)
        a.set_title(lab, loc="left"); a.set_xlabel("distance from the separatrix, mapped to the midplane [cm]")
    axes[0, 0].legend(frameon=False, fontsize=8.5)
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)

profiles([("FHRP", "density", "Midplane $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("FHRP", "electron_temp", "Midplane $T_e$ [eV]", 1), ("FHRP", "jsat", "Midplane $J_{sat}$ [kA m$^{-2}$]", 1e-3),
          ("TS", "density", "Divertor entrance $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("TS", "electron_temp", "Divertor entrance $T_e$ [eV]", 1)],
         out + "/upstream.png", 3, "Upstream: outboard midplane probe and divertor-entrance Thomson scattering")
profiles([("LFS-LP", "density", "Outer target $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("LFS-LP", "electron_temp", "Outer target $T_e$ [eV]", 1),
          ("LFS-LP", "jsat", "Outer target $J_{sat}$ [kA m$^{-2}$]", 1e-3), ("LFS-IR", "q_parallel", "Outer target $q_\\parallel$ [MW m$^{-2}$]", 1e-6)],
         out + "/outer_target.png", 2, "Outer target (floor): Langmuir probes and infrared camera")
profiles([("HFS-LP", "density", "Inner target $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("HFS-LP", "electron_temp", "Inner target $T_e$ [eV]", 1),
          ("HFS-LP", "jsat", "Inner target $J_{sat}$ [kA m$^{-2}$]", 1e-3)],
         out + "/inner_target.png", 3, "Inner target (inner wall): Langmuir probes")
print("done")
