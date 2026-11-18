"""Balances and regional ionisation of a Hermes-3 run on the TCV-X21 GOAT grid.

usage (from us_runs/): python tcvx21/analyse.py <run dir> [fields.png]
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata
import compare as c

run = sys.argv[1]
r = usdata.Run(run, c.GRID); g = r.g; n = r.n; it = -1; twopi = 2*np.pi; qe = c.qe
Ni = r("Nd+", True); Vi = r("Vd+", True); Te = r("Te", True); Ti = r("Td+", True); Nn = r("Nd", True); Pe = r("Pe", True); Pi = r("Pd+", True)
L, R, kind = g["face_L"], g["face_R"], g["face_kind"]; fb = g["face_flux"]/g["face_B"]; wall, core = kind == 2, kind == 1
psin = c.psin_RZ(g["R"], g["Z"])
print("t = %.1f ms, %d outputs" % (1e3*r.t[it], len(r.t)))
cont = r.volume_integral(Ni); cn = r.volume_integral(Nn); W = r.volume_integral(1.5*(Pe + Pi))
k = min(5, len(r.t) - 1); dt = r.t[-1] - r.t[-1-k]
print("ion content %.3e, neutral content %.3e, stored energy %.2f kJ | rates of change over the last %.0f ms: ions %.2e /s, neutrals %.2e /s, relative %.1e and %.1e per 10 ms" % (
    cont[-1], cn[-1], W[-1]/1e3, 1e3*dt, (cont[-1]-cont[-1-k])/dt, (cn[-1]-cn[-1-k])/dt, (cont[-1]/cont[-1-k]-1)*0.01/dt, (cn[-1]/cn[-1-k]-1)*0.01/dt))
sign = np.where(fb >= 0, 1.0, -1.0)
ion_wall = twopi*(Ni[it][R]*np.maximum(sign*Vi[it][R], 0)*np.abs(fb))[wall].sum()
D = 0.2
core_in = twopi*(D*g["face_tiso"]*(Ni[it][R] - Ni[it][L]))[core].sum()
out = "particles [1/s]: ion flow to walls %.3e | diffusive ion inflow through the core boundary %.3e" % (ion_wall, core_in)
for name, lab in (("Sd_puff", "puff"), ("Sd_pump", "neutrals absorbed by the wall"), ("Sd_target_recycle", "recycled neutrals"), ("Sd+_src", "ion source"), ("Sd_src", "neutral source")):
    if r.has(name):
        out += " | %s %.3e" % (lab, r.volume_integral(r(name, True)[it]))
print(out)
P_in = sum(1.5*r.volume_integral(r("P%s_src" % s, True)[it]) for s in ("d+", "e"))
P_sh = -(r.volume_integral(r("Ee_sheath", True)[it]) + r.volume_integral(r("Ei_sheath", True)[it]))
print("power [kW]: heating %.1f | sheath %.1f | 13.6 eV per ion to the wall %.1f | remainder %.1f" % (P_in/1e3, P_sh/1e3, 13.6*qe*ion_wall/1e3, (P_in - P_sh)/1e3))
# regions as in Wang et al.: core (inside the separatrix above the X-point), LFS and HFS divertor, rest of the SOL
Rc, Zc = g["R"][:n], g["Z"][:n]
ix = np.argmin(np.where((np.abs(psin[:n] - 1) < 0.01) & (Zc < -0.2), np.hypot(g["Bp"][:n], 0), np.inf)); Rx, Zx = Rc[ix], Zc[ix]
closed = (psin[:n] < 1) & (Zc > Zx)
div = Zc < Zx
regions = {"inside the separatrix": closed, "LFS divertor (below the X-point, R > R_X)": div & (Rc > Rx), "HFS divertor (below the X-point, R < R_X)": div & (Rc <= Rx), "rest of the SOL": ~closed & ~div}
print("X-point (lowest poloidal field near the separatrix) at R %.3f, Z %.3f" % (Rx, Zx))
izname = [v for v in r.d.variables if v.startswith("Sd+_iz") or v == "Sd+_iz"]
names = [v for v in r.d.variables if r.d[v].ndim == 4]
iz = [v for v in names if "iz" in v.lower() and v.startswith("S")]
print("reaction diagnostics present:", [v for v in names if v[0] in "SRE" and any(k in v for k in ("iz", "rec", "cx", "ex"))])
if iz:
    S = np.abs(r(iz[0])[it]); tot = twopi*(S*g["vol"][:n]).sum()
    print("ionisation %.3e /s (%s): " % (tot, iz[0]) + ", ".join("%s %.1f%%" % (k, 100*twopi*(S*g["vol"][:n])[m].sum()/tot) for k, m in regions.items()))
    print("   SOLPS-ITER (Wang et al.): 65.0% in the SOL, of which LFS divertor 31.2% and HFS divertor 14.9% of the total")
rad = [v for v in names if v.startswith("R") and not v.startswith("resid")]
for v in rad:
    print("radiation %s: %.1f kW" % (v, abs(r.volume_integral(r(v, True)[it]))/1e3))
omp = np.flatnonzero((np.abs(Zc) < 0.012) & (Rc > c.q["rmagx"])); omp = omp[np.argsort(psin[omp])]
x = c.rsep_of_psin(psin[omp])
f = lambda a: float(np.interp(0.0, x, a[it][omp]))
sol = omp[(x[np.arange(len(omp))] > 0) & (x < 1.5)]
lam = lambda a: -10/np.polyfit(c.rsep_of_psin(psin[sol]), np.log(a[it][sol]), 1)[0]
print("midplane separatrix: n_e %.2e m^-3, Te %.1f eV, Ti %.1f eV | decay lengths over the first 1.5 cm: n %.1f mm, Te %.1f mm | first cell (%.2f cm): n %.2e, Te %.0f, Ti %.0f" % (
    f(Ni), f(Te), f(Ti), lam(Ni), lam(Te), x[0], Ni[it][omp[0]], Te[it][omp[0]], Ti[it][omp[0]]))
print("neutral density: max %.2e, median %.2e, at the midplane wall cell %.2e m^-3" % (Nn[it][:n].max(), np.median(Nn[it][:n]), Nn[it][omp[-1]]))
if len(sys.argv) > 2:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots(1, 4, figsize=(17, 8))
    panels = [(Ni[it, :n], "electron density [m$^{-3}$]", dict(norm=LogNorm(1e16, Ni[it, :n].max()), cmap="viridis")),
              (Te[it, :n], "electron temperature [eV]", dict(norm=LogNorm(0.5, Te[it, :n].max()), cmap="inferno")),
              (Nn[it, :n], "D neutral density [m$^{-3}$]", dict(norm=LogNorm(1e14, Nn[it, :n].max()), cmap="viridis"))]
    if iz:
        S = np.abs(r(iz[0])[it]); panels.append((np.maximum(S, 1.0), "ionisation source [m$^{-3}$ s$^{-1}$]", dict(norm=LogNorm(S.max()*1e-4, S.max()), cmap="magma")))
    else:
        panels.append((Ti[it, :n], "ion temperature [eV]", dict(norm=LogNorm(0.5, Ti[it, :n].max()), cmap="inferno")))
    for a, (val, lab, kw) in zip(ax, panels):
        pc = usdata.plot_cells(a, g, val, edgecolors="none", **kw); fig.colorbar(pc, ax=a, shrink=0.8); a.set_title(lab); a.set_xlabel("R [m]")
    ax[0].set_ylabel("Z [m]")
    fig.suptitle(sys.argv[3] if len(sys.argv) > 3 else run)
    fig.savefig(sys.argv[2], dpi=100, bbox_inches="tight")
