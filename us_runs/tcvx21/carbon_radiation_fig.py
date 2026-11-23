"""Carbon radiation in the poloidal plane for several carbon runs, on one colour scale.

usage (from us_runs/): python tcvx21/carbon_radiation_fig.py <out.png> <label>=<run dir> [...]

The first run is the reference: the others are also shown as a ratio to it.
The radiation is rebuilt from the output with the OpenADAS tables the run uses
(line radiation of C..C3+, recombination and bremsstrahlung onto C+..C4+).
"""
import sys, json
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata
import compare as c
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

INK, MUTED = "#0b0b0b", "#898781"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e", "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.titlesize": 10.5, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
WALL = np.array([(0.671, 0.750), (0.624, 0.700), (0.624, -0.700), (0.671, -0.750), (0.965, -0.750), (1.136, -0.543), (1.136, 0.543), (0.965, 0.750), (0.671, 0.750)])
twopi = 2*np.pi


def adas(kind):
    d = json.load(open("json_database/%s96_c.json" % kind))
    lt, ln, co = np.array(d["log_temperature"]), np.array(d["log_density"]), np.array(d["log_coeff"])

    def f(level, te, n_e):
        x = np.clip(np.log10(np.maximum(te, 1e-30)), lt[0], lt[-1]); y = np.clip(np.log10(np.maximum(n_e, 1.0)), ln[0], ln[-1])
        i = np.clip(np.searchsorted(lt, x), 1, len(lt) - 1); j = np.clip(np.searchsorted(ln, y), 1, len(ln) - 1)
        wx = (x - lt[i - 1])/(lt[i] - lt[i - 1]); wy = (y - ln[j - 1])/(ln[j] - ln[j - 1])
        t = co[level]
        return 10**(t[i - 1, j - 1]*(1 - wx)*(1 - wy) + t[i, j - 1]*wx*(1 - wy) + t[i - 1, j]*(1 - wx)*wy + t[i, j]*wx*wy)
    return f


plt_c, prb_c = adas("plt"), adas("prb")
out = sys.argv[1]
runs = [a.split("=", 1) for a in sys.argv[2:]]
data = []
for label, run in runs:
    r = usdata.Run(run, c.GRID); g = r.g; n = r.n; vol = g["vol"][:n]
    ions = [s for s in ["c+"] + ["c+%d" % k for k in range(2, 7)] if r.has("N" + s)]
    Te = r("Te")[-1]; N = {s: r("N" + s)[-1] for s in ["c"] + ions}
    ne = r("Nd+")[-1] + sum((k + 1)*N[s] for k, s in enumerate(ions))
    states = ["c"] + ions
    by_state = {}
    prad = np.zeros(n)
    for k in range(len(ions)):
        line = ne*N[states[k]]*plt_c(k, Te, ne); rec = ne*N[states[k + 1]]*prb_c(k, Te, ne)
        by_state[states[k]] = by_state.get(states[k], 0) + twopi*(line*vol).sum()
        by_state[states[k + 1]] = by_state.get(states[k + 1], 0) + twopi*(rec*vol).sum()
        prad += line + rec
    psin = c.psin_RZ(g["R"][:n], g["Z"][:n]); Z = g["Z"][:n]; R = g["R"][:n]
    Zx = -0.41
    regions = {"inside the separatrix": (psin < 1) & (Z > Zx), "outer divertor (below the X-point, R > 0.75)": (Z < Zx) & (R > 0.75),
               "inner divertor (below the X-point, R < 0.75)": (Z < Zx) & (R <= 0.75), "rest of the SOL": ~((psin < 1) & (Z > Zx)) & (Z >= Zx)}
    tot = twopi*(prad*vol).sum()
    print("%s: %.2f kW | " % (label, tot/1e3) + ", ".join("%s %.2f kW" % (k, twopi*(prad*vol)[m].sum()/1e3) for k, m in regions.items()))
    print("     by state [kW]: " + ", ".join("%s %.2f" % (k, v/1e3) for k, v in by_state.items()) + " | peak %.2e W/m3" % prad.max())
    data.append((label, prad, tot))
polys = [np.stack([g["vertex_R"][v[v >= 0]], g["vertex_Z"][v[v >= 0]]], 1) for v in g["cell_vertices"]][:n]
RR, ZZ = np.meshgrid(np.linspace(c.er[0], c.er[-1], 300), np.linspace(c.ez[0], c.ez[-1], 600), indexing="ij")
PSN = c.psin_RZ(RR, ZZ)

nref = len(data); nratio = len(data) - 1
fig, ax = plt.subplots(1, nref + nratio, figsize=(3.05*(nref + nratio), 7.4))
vmax = max(d[1].max() for d in data); norm = LogNorm(vmax*1e-3, vmax)
for a, (label, prad, tot) in zip(ax[:nref], data):
    pc = PolyCollection(polys, array=np.clip(prad, norm.vmin, None), norm=norm, cmap="Blues", edgecolors="face", linewidths=0.1)
    a.add_collection(pc); a.set_title("%s\n%.1f kW" % (label, tot/1e3))
cb = fig.colorbar(pc, ax=list(ax[:nref]), orientation="horizontal", fraction=0.04, pad=0.08, aspect=40); cb.set_label("carbon radiation [W m$^{-3}$]"); cb.outline.set_visible(False)
for a, (label, prad, tot) in zip(ax[nref:], data[1:]):
    ratio = prad/np.maximum(data[0][1], 1e-30)
    shown = data[0][1] > norm.vmin                       # only where there is radiation to compare
    pc2 = PolyCollection([p for p, s in zip(polys, shown) if s], array=np.log2(ratio[shown]), cmap="RdBu_r", clim=(-1, 1), edgecolors="face", linewidths=0.1)
    a.add_collection(PolyCollection([p for p, s in zip(polys, shown) if not s], facecolors="#f0efec", edgecolors="none"))
    a.add_collection(pc2); a.set_title("%s\n/ %s" % (label, data[0][0]))
cb2 = fig.colorbar(pc2, ax=list(ax[nref:]), orientation="horizontal", fraction=0.04, pad=0.08, aspect=27, ticks=[-1, -0.585, 0, 0.585, 1])
cb2.ax.set_xticklabels(["0.5", "0.67", "1", "1.5", "2"]); cb2.set_label("ratio of carbon radiation"); cb2.outline.set_visible(False)
for k, a in enumerate(ax):
    a.plot(WALL[:, 0], WALL[:, 1], color=INK, lw=1.2, zorder=5)
    a.contour(RR, ZZ, PSN, levels=[1.0], colors=[INK], linewidths=0.7, linestyles="--", zorder=4)
    a.set_aspect("equal"); a.set_xlim(0.59, 1.16); a.set_ylim(-0.78, 0.78); a.set_xlabel("R [m]")
    if k: a.set_yticklabels([])
ax[0].set_ylabel("Z [m]")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)
