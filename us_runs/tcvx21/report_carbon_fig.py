"""Carbon figure for the report: carbon ion density, Z_eff and carbon radiation.

usage (from us_runs/): python tcvx21/report_carbon_fig.py <carbon run dir> <out.png>
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import json
import numpy as np, usdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, Normalize
import compare as c

run, out = sys.argv[1:3]
r = usdata.Run(run, c.GRID); g = r.g; n = r.n
ions = [s for s in ["c+"] + ["c+%d" % k for k in range(2, 7)] if r.has("N" + s)]
Nc = {s: r("N" + s)[-1] for s in ["c"] + ions}
Nd = r("Nd+")[-1]; Te = r("Te")[-1]
nz = sum(Nc[s] for s in ions)
ne = Nd + sum((k + 1)*Nc[s] for k, s in enumerate(ions))
zeff = (Nd + sum((k + 1)**2*Nc[s] for k, s in enumerate(ions)))/ne

def table(kind):
    d = json.load(open("json_database/%s96_c.json" % kind))
    lt, ln, co = np.array(d["log_temperature"]), np.array(d["log_density"]), np.array(d["log_coeff"])
    def f(level):
        x = np.clip(np.log10(np.maximum(Te, 1e-30)), lt[0], lt[-1]); y = np.clip(np.log10(np.maximum(ne, 1.0)), ln[0], ln[-1])
        i = np.clip(np.searchsorted(lt, x), 1, len(lt) - 1); j = np.clip(np.searchsorted(ln, y), 1, len(ln) - 1)
        wx = (x - lt[i - 1])/(lt[i] - lt[i - 1]); wy = (y - ln[j - 1])/(ln[j] - ln[j - 1]); a = co[level]
        return 10**(a[i - 1, j - 1]*(1 - wx)*(1 - wy) + a[i, j - 1]*wx*(1 - wy) + a[i - 1, j]*(1 - wx)*wy + a[i, j]*wx*wy)
    return f
plt_c, prb_c = table("plt"), table("prb"); states = ["c"] + ions
prad = sum(ne*Nc[states[k]]*plt_c(k) + ne*Nc[states[k + 1]]*prb_c(k) for k in range(len(ions)))

polys = [np.stack([g["vertex_R"][v[v >= 0]], g["vertex_Z"][v[v >= 0]]], 1) for v in g["cell_vertices"]][:n]
WALL = np.array([(0.671, 0.750), (0.624, 0.700), (0.624, -0.700), (0.671, -0.750), (0.965, -0.750), (1.136, -0.543), (1.136, 0.543), (0.965, 0.750), (0.671, 0.750)])
RR, ZZ = np.meshgrid(np.linspace(c.er[0], c.er[-1], 300), np.linspace(c.ez[0], c.ez[-1], 600), indexing="ij")
PSN = c.psin_RZ(RR, ZZ)
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e", "xtick.color": "#898781", "ytick.color": "#898781",
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 3, figsize=(10.5, 6.6))
zmask = np.where(ne > 3e17, zeff, np.nan)       # Z_eff means little where there is almost no plasma
panels = [(np.maximum(nz, 1e14), "carbon ion density [m$^{-3}$]", LogNorm(1e15, 1e18), "viridis"),
          (zmask, "Z$_{eff}$ (where $n_e$ > 3$\\times$10$^{17}$ m$^{-3}$)", Normalize(1.0, 2.5), "viridis"),
          (np.maximum(prad, 1.0), "carbon radiation [W m$^{-3}$]", LogNorm(1e3, 1e6), "inferno")]
for a, (val, lab, norm, cmap) in zip(ax, panels):
    cm = plt.get_cmap(cmap).copy(); cm.set_bad("#f0efec")
    pc = PolyCollection(polys, array=np.ma.masked_invalid(val), norm=norm, cmap=cm, edgecolors="face", linewidths=0.1)
    a.add_collection(pc)
    a.plot(WALL[:, 0], WALL[:, 1], color="#0b0b0b", lw=1.3)
    a.contour(RR, ZZ, PSN, levels=[1.0], colors=["#0b0b0b"], linewidths=0.7, linestyles="--")
    a.set_aspect("equal"); a.set_xlim(0.59, 1.16); a.set_ylim(-0.78, 0.78); a.set_xlabel("R [m]")
    cb = fig.colorbar(pc, ax=a, orientation="horizontal", fraction=0.045, pad=0.09); cb.set_label(lab); cb.outline.set_visible(False)
ax[0].set_ylabel("Z [m]")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)
