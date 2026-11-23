"""Potential, currents and drift flows of a TCV-X21 run with drifts_unstructured.

usage (from us_runs/): python tcvx21/drifts.py <run dir> [<reference run dir>] [figure.png]

Prints the potential at the midplane and the targets, the current balance,
the radial electric field at the midplane against two estimates, and, if a
reference run is given, how far density and temperatures are from it.
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata
import compare as c

args = [a for a in sys.argv[1:] if not a.endswith(".png")]
fig = [a for a in sys.argv[1:] if a.endswith(".png")]
run = args[0]
grid = c.GRID
if "@@" in run:
    run, grid = run.split("@@")
r = usdata.Run(run, grid); g = r.g; n = r.n; it = -1; twopi = 2*np.pi; qe = c.qe
Ni = r("Nd+", True)[it]; Vi = r("Vd+", True)[it]; Te = r("Te", True)[it]; Ti = r("Td+", True)[it]; Pi = r("Pd+", True)[it]
phi = r("phi", True)[it]
L, R, kind = g["face_L"], g["face_R"], g["face_kind"]; fb = g["face_flux"]/g["face_B"]; wall = kind == 2
psin = c.psin_RZ(g["R"], g["Z"]); Rc, Zc = g["R"][:n], g["Z"][:n]
print("%s: t = %.2f ms" % (run, 1e3*r.t[it]))

# ---- steadiness
if len(r.t) > 2:
    k = min(5, len(r.t) - 1); dt = r.t[-1] - r.t[-1-k]
    for name in ("Nd+", "Nd", "Pe", "phi"):
        a = r(name)
        w = g["vol"][:n]
        change = np.sqrt((w*(a[-1] - a[-1-k])**2).sum()/(w*a[-1]**2).sum())
        print("   %-4s relative change over the last %.1f ms: %.2e" % (name, 1e3*dt, change))

# ---- midplane profile
omp = np.flatnonzero((np.abs(Zc) < 0.012) & (Rc > c.q["rmagx"])); omp = omp[np.argsort(psin[omp])]
x = c.rsep_of_psin(psin[omp])                         # cm
sep = lambda a: float(np.interp(0.0, x, a[omp]))
print("midplane separatrix: n_e %.2e m^-3, Te %.1f eV, Ti %.1f eV, phi %.1f V (%.2f Te)" % (sep(Ni), sep(Te), sep(Ti), sep(phi), sep(phi)/sep(Te)))
Rm = Rc[omp]
Er = -np.gradient(phi[omp], Rm)
Er_sol = 3*(-np.gradient(Te[omp], Rm))                # -3 dTe/dr, the sheath-connected estimate
Er_fb = np.gradient(Pi[omp], Rm)/(qe*Ni[omp])         # ion force balance without flow
print("midplane radial electric field [kV/m]  (x = distance from the separatrix)")
print("   x [cm]    E_r   -3 dTe/dr   dpi/dr/(e n)    phi [V]   Te [eV]")
for j in range(0, len(omp), max(1, len(omp)//14)):
    print("   %6.2f  %6.2f   %6.2f      %6.2f      %7.1f   %6.1f" % (x[j], Er[j]/1e3, Er_sol[j]/1e3, Er_fb[j]/1e3, phi[omp[j]], Te[omp[j]]))

# ---- targets: sheath-edge values in the ghost cells
sign = np.where(fb >= 0, 1.0, -1.0)
Ve = r("Ve", True)[it]
gam_i = np.where(wall, Ni[R]*np.maximum(sign*Vi[R], 0), 0.0)
gam_e = np.where(wall, Ni[R]*np.maximum(sign*Ve[R], 0), 0.0)
jsh = qe*(gam_i - gam_e)                              # current density along B into the wall
I_face = twopi*jsh*np.abs(fb)
fR, fZ = g["face_Rmid"], g["face_Zmid"]
regions = {"outer target (floor, R > 0.69)": wall & (fZ < -0.749) & (fR > 0.69),
           "inner target (inner wall, Z < -0.2)": wall & (fR < 0.6245) & (fZ < -0.2),
           "rest of the wall": wall & ~((fZ < -0.749) & (fR > 0.69)) & ~((fR < 0.6245) & (fZ < -0.2))}
print("sheath currents (positive = ions in excess, current into the wall)")
for name, m in regions.items():
    sat = twopi*(qe*gam_i*np.abs(fb))[m].sum()
    peak = np.argmax(np.where(m, gam_i, 0))
    print("   %-36s net %8.2f A | ion saturation current %8.1f A | at the peak ion flux: Te %.1f eV, phi %.1f V (%.2f Te), j/j_sat %+.2f" % (
        name, I_face[m].sum(), sat, Te[R[peak]], phi[L[peak]], phi[L[peak]]/max(Te[R[peak]], 1e-3), jsh[peak]/max(qe*gam_i[peak], 1e-30)))
print("   total current into the wall %.3f A (should vanish in steady state; ion saturation current %.1f A)" % (I_face[wall].sum(), twopi*(qe*gam_i*np.abs(fb))[wall].sum()))
if r.has("DivJpar"):
    for name in ("DivJpar", "DivJsheath", "DivJdrift", "DivJperp"):
        a = r(name)[it]
        print("   integral of |%s| %.2f A, net %.3f A" % (name, twopi*(np.abs(a)*g["vol"][:n]).sum(), twopi*(a*g["vol"][:n]).sum()))
    tot = sum(r(name)[it] for name in ("DivJpar", "DivJsheath", "DivJdrift", "DivJperp"))
    print("   integral of |Div J| %.3f A" % (twopi*(np.abs(tot)*g["vol"][:n]).sum()))
if r.has("Se_exb"):
    for name, lab in (("Se_exb", "E x B"), ("Se_magdrift", "magnetic drift of the electrons")):
        a = r(name)[it]
        print("   particles moved by the %s: %.2e /s (integral of |source|), net %.2e /s" % (lab, twopi*(np.abs(a)*g["vol"][:n]).sum(), twopi*(a*g["vol"][:n]).sum()))
ion_wall = twopi*(gam_i*np.abs(fb))[wall].sum()
print("ion flow to the walls %.3e /s: outer target %.3e, inner target %.3e" % (ion_wall, twopi*(gam_i*np.abs(fb))[regions["outer target (floor, R > 0.69)"]].sum(),
      twopi*(gam_i*np.abs(fb))[regions["inner target (inner wall, Z < -0.2)"]].sum()))

# ---- against a reference run
if len(args) > 1:
    ref_dir, ref_grid = (args[1].split("@@") + [grid])[:2]
    ref = usdata.Run(ref_dir, ref_grid)
    print("against %s (t = %.1f ms), real cells with n_e > 1e17 m^-3:" % (ref_dir, 1e3*ref.t[-1]))
    if ref.n != n:
        print("   (different grids: skipped)")
        sys.exit(0)
    m = Ni[:n] > 1e17
    for name in ("Nd+", "Te", "Td+", "Nd"):
        a, b = r(name)[it], ref(name)[-1]
        rel = np.abs(a/b - 1)[m]
        print("   %-4s median %.2f%%, 95th percentile %.2f%%, max %.2f%%" % (name, 100*np.median(rel), 100*np.percentile(rel, 95), 100*rel.max()))

if fig:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 3, figsize=(13, 6.5))
    Jp = r("Jpar")[it] if r.has("Jpar") else np.zeros(n)
    for a, val, lab, kw in ((ax[0], phi[:n], "potential [V]", dict(cmap="viridis")),
                            (ax[1], phi[:n]/np.maximum(Te[:n], 0.5), "potential / T_e", dict(cmap="viridis", clim=(0, 5))),
                            (ax[2], Jp/1e3, "parallel current density [kA m$^{-2}$]", dict(cmap="RdBu_r", clim=(-np.percentile(np.abs(Jp), 98)/1e3, np.percentile(np.abs(Jp), 98)/1e3)))):
        pc = usdata.plot_cells(a, g, val, **kw); f.colorbar(pc, ax=a, orientation="horizontal", pad=0.07, label=lab); a.set_xlabel("R [m]")
    ax[0].set_ylabel("Z [m]")
    f.savefig(fig[0], dpi=130, bbox_inches="tight")
    print("wrote", fig[0])
