"""Carbon in a TCV-X21 GOAT-grid run with wall sputtering.

usage (from us_runs/): python tcvx21/carbon.py <run dir> [figure.png]
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata
import compare as c

run = sys.argv[1]
r = usdata.Run(run, c.GRID); g = r.g; n = r.n; it = -1; twopi = 2*np.pi; qe = c.qe
ions = [s for s in ["c+"] + ["c+%d" % k for k in range(2, 7)] if r.has("N" + s)]   # charge states in the run
ne = r("Ne")[it] if r.has("Ne") else None
Nd = r("Nd+")[it]
Nc = {s: r("N" + s)[it] for s in ["c"] + ions}
nz = sum(Nc[s] for s in ions)
if ne is None:
    ne = Nd + sum((k + 1)*Nc[s] for k, s in enumerate(ions))
zeff = (Nd + sum((k + 1)**2*Nc[s] for k, s in enumerate(ions)))/ne
vol = g["vol"][:n]
integ = lambda a: twopi*(a*vol).sum()
print("t = %.1f ms, %d outputs" % (1e3*r.t[it], len(r.t)))
tot = sum(integ(Nc[s]) for s in ions)
hist = [sum(twopi*(r("N" + s)[k]*vol).sum() for s in ions) for k in range(len(r.t))]
if len(hist) > 2:
    print("carbon ion content %.3e; change over the last output %.1f%%" % (tot, 100*(hist[-1]/hist[-2] - 1)))
print("content: atoms %.3e | " % integ(Nc["c"]) + " | ".join("%s %.2e" % (s, integ(Nc[s])) for s in ions))
print("carbon ions / electrons, whole volume: %.2f%%" % (100*tot/integ(ne)))
if r.has("Sc_sputter"):
    print("sputtered atoms: %.3e /s" % integ(r("Sc_sputter")[it]))
for name, lab in (("Sc_pump", "carbon atoms absorbed by the wall"),):
    if r.has(name):
        print("%s: %.3e /s" % (lab, -integ(r(name)[it])))
psin = c.psin_RZ(g["R"][:n], g["Z"][:n])
omp = np.flatnonzero((np.abs(g["Z"][:n]) < 0.012) & (g["R"][:n] > c.q["rmagx"])); omp = omp[np.argsort(psin[omp])]
x = c.rsep_of_psin(psin[omp])
sep = lambda a: float(np.interp(0.0, x, a[omp]))
print("midplane separatrix: n_e %.2e, Te %.1f eV, carbon concentration %.2f%%, Zeff %.2f | core edge (first cell): concentration %.2f%%, Zeff %.2f" % (
    sep(ne), sep(r("Te")[it]), 100*sep(nz/ne), sep(zeff), 100*(nz/ne)[omp[0]], zeff[omp[0]]))


def adas(kind):
    """OpenADAS loss coefficient [W m^3] of each charge state as a function (Te [eV], ne [m^-3]),
    interpolated as in src/adas_reaction.cxx (bilinear in log10, clamped to the table)"""
    import json
    d = json.load(open("json_database/%s96_c.json" % kind))
    lt, ln, co = np.array(d["log_temperature"]), np.array(d["log_density"]), np.array(d["log_coeff"])

    def f(level, te, n_e):
        x = np.clip(np.log10(np.maximum(te, 1e-30)), lt[0], lt[-1]); y = np.clip(np.log10(np.maximum(n_e, 1.0)), ln[0], ln[-1])
        i = np.clip(np.searchsorted(lt, x), 1, len(lt) - 1); j = np.clip(np.searchsorted(ln, y), 1, len(ln) - 1)
        wx = (x - lt[i - 1])/(lt[i] - lt[i - 1]); wy = (y - ln[j - 1])/(ln[j] - ln[j - 1])
        c_ = co[level]
        return 10**(c_[i - 1, j - 1]*(1 - wx)*(1 - wy) + c_[i, j - 1]*wx*(1 - wy) + c_[i - 1, j]*(1 - wx)*wy + c_[i, j]*wx*wy)
    return f


# Carbon radiation, rebuilt from the output with the same OpenADAS tables the run uses:
# line radiation of each state (plt) and recombination plus bremsstrahlung (prb)
Te = r("Te")[it]
plt_c, prb_c = adas("plt"), adas("prb")
states = ["c"] + ions
prad = np.zeros(n)
for k in range(len(ions)):
    prad += ne*Nc[states[k]]*plt_c(k, Te, ne) + ne*Nc[states[k + 1]]*prb_c(k, Te, ne)
# (line radiation of the highest state followed is not in the run: it has no ionisation reaction)
total_c = integ(prad)
closed = (psin < 1) & (g["Z"][:n] > -0.42)
print("carbon radiation (line + recombination + bremsstrahlung): %.1f kW = %.0f%% of 160 kW; inside the separatrix %.1f kW, below Z = -0.42 m %.1f kW" % (
    total_c/1e3, 100*total_c/160e3, integ(prad*closed)/1e3, integ(prad*(g["Z"][:n] < -0.42))/1e3))
for v in [v for v in r.d.variables if r.d[v].ndim == 4 and v.startswith("R")]:
    print("radiation diagnostic %-10s %.2f kW" % (v, abs(integ(r(v)[it]))/1e3))
rad = []
if len(sys.argv) > 2:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots(1, 4, figsize=(17, 8))
    panels = [(np.maximum(nz, 1e10), "carbon ion density [m$^{-3}$]", LogNorm(1e15, max(nz.max(), 2e15)), "viridis"),
              (np.maximum(Nc["c"], 1e10), "carbon atom density [m$^{-3}$]", LogNorm(1e14, max(Nc["c"].max(), 2e14)), "viridis"),
              (zeff, "Z$_{eff}$", None, "magma"),
              (np.maximum(prad, 1.0), "carbon radiation [W m$^{-3}$]", LogNorm(max(prad.max(), 10)*1e-3, max(prad.max(), 10)), "inferno")]
    for a, (val, lab, norm, cmap) in zip(ax, panels):
        pc = usdata.plot_cells(a, g, val, edgecolors="none", norm=norm, cmap=cmap); fig.colorbar(pc, ax=a, shrink=0.8); a.set_title(lab); a.set_xlabel("R [m]")
    ax[0].set_ylabel("Z [m]")
    fig.savefig(sys.argv[2], dpi=100, bbox_inches="tight")
