"""Compare a Hermes-3 run on the TCV-X21 GOAT grid with the TCV-X21 experimental dataset.

usage: python tcvx21/compare.py <figure.png> <label>=<run dir> [<label>=<run dir> ...]
       (run from us_runs/)
       A run on a structured (hypnotoad) grid is given as <label>=<run dir>::<grid file>.

Observables (all against R^u - R^u_sep, the distance from the separatrix mapped
along flux surfaces to the outboard midplane, as in the dataset):
  FHRP    outboard midplane (Z = 0): n_e, T_e, J_sat
  TS      divertor entrance (Thomson chord at R = 0.9 m, below the plasma): n_e, T_e
  LFS-LP  floor target (Z = -0.75 m): n_e, T_e, J_sat
  HFS-LP  inner-wall target (R = 0.624 m): n_e, T_e, J_sat
  LFS-IR  floor target: parallel heat flux

For each observable the normalised distance of the dataset is printed,
    d_j = sqrt( mean_i (s_i - e_i)^2 / de_i^2 ),
with the simulation interpolated linearly onto the experimental positions
inside its own range, next to the values of the SOLPS-ITER "standard" case
(Wang et al., Nucl. Fusion 64 (2024) 056040, table II).

Target values are the sheath-edge values that the sheath boundary stores in
the ghost cell of each wall face. J_sat = e * (parallel ion flux density).
q_parallel = sheath heat flux + 13.6 eV per ion, along the field.
"""
import sys
sys.path.insert(0, "../tools/unstructured")
import numpy as np, netCDF4, usdata
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq
from hypnotoad.geqdsk._geqdsk import read as read_geqdsk
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

GRID = "grids/tcvx21_goat.grd.nc"
DATA = "tcvx21/data/"
qe = 1.602176634e-19; twopi = 2*np.pi
GAMMA_E, GAMMA_I = 4.5, 2.5

# SOLPS-ITER standard case, d_j against forward / reversed field data (table II of Wang et al.)
SOLPS = {("FHRP", "density"): (0.800, 1.275), ("FHRP", "electron_temp"): (0.359, 0.907), ("FHRP", "jsat"): (1.523, 1.845),
         ("TS", "density"): (0.848, 1.148), ("TS", "electron_temp"): (0.818, 1.095),
         ("LFS-LP", "density"): (3.222, 2.046), ("LFS-LP", "electron_temp"): (2.170, 1.907), ("LFS-LP", "jsat"): (3.688, 2.342),
         ("HFS-LP", "density"): (3.360, 4.269), ("HFS-LP", "electron_temp"): (3.136, 2.750), ("HFS-LP", "jsat"): (2.079, 3.702),
         ("LFS-IR", "q_parallel"): (6.350, 6.322)}

# ---- equilibrium: psiN(R, Z) and the midplane mapping
with open(DATA + "1.experimental_data_reference_scenario_65402_t1.eqdsk") as f:
    q = read_geqdsk(f)
er = q["rleft"] + q["rdim"]*np.arange(q["nx"])/(q["nx"] - 1)
ez = q["zmid"] - q["zdim"]/2 + q["zdim"]*np.arange(q["ny"])/(q["ny"] - 1)
spl = RectBivariateSpline(er, ez, (np.asarray(q["psi"]) - q["simagx"])/(q["sibdry"] - q["simagx"]), kx=3, ky=3)
psin_RZ = lambda R, Z: spl(R, Z, grid=False)
Z_OMP = 0.0
R_sep = brentq(lambda R: psin_RZ(R, Z_OMP) - 1.0, q["rmagx"] + 0.05, er[-1] - 1e-3)
_Rm = np.linspace(q["rmagx"] + 0.02, er[-1], 4000); _pm = psin_RZ(_Rm, Z_OMP)
assert np.all(np.diff(_pm) > 0)
rsep_of_psin = lambda p: 100*(np.interp(p, _pm, _Rm) - R_sep)      # cm
psin_of_rsep = lambda r: np.interp(R_sep + np.asarray(r)/100, _Rm, _pm)


def load(run):
    r = usdata.Run(run, GRID); g = r.g; n = r.n; it = -1
    Ni = r("Nd+", True)[it]; Vi = r("Vd+", True)[it]; Ti = r("Td+", True)[it]; Te = r("Te", True)[it]
    L, R = g["face_L"], g["face_R"]; wall = g["face_kind"] == 2
    fb = g["face_flux"]/g["face_B"]; sign = np.where(fb >= 0, 1.0, -1.0)
    gam = np.where(wall, Ni[R]*np.maximum(sign*Vi[R], 0), 0.0)          # parallel ion flux density at the sheath edge
    qpar = np.where(wall, (GAMMA_E*Te[R] + GAMMA_I*Ti[R])*qe*gam, 0.0)
    psin = psin_RZ(g["R"], g["Z"])                                       # cells and ghosts (ghosts sit on their face)
    grid_psin = 1 + g["psi"]*twopi/q["simagx"]   # grid file: Wb/rad, sign flipped; this eqdsk: Wb
    return dict(r=r, g=g, n=n, t=r.t[it], Ne=Ni, Te=Te, Ti=Ti, gam=gam, qpar=qpar, psin=psin, grid_psin=grid_psin, wall=wall, L=L, R=R, fb=fb)


def load_structured(run, gridfile):
    """Profiles of a run on a structured grid (NXPE = 1), for the same observables.

    Target values are those on the target boundary: the mean of the last cell
    and its guard cell. Each x index is a flux surface, so its psiN is taken
    from the last row of the outer leg, which lies inside the geqdsk box.
    """
    import glob
    files = sorted(glob.glob(run + "/BOUT.dmp.*.nc"), key=lambda f: int(f.split(".")[-2]))
    ds = [netCDF4.Dataset(f) for f in files]
    myg = int(ds[0]["MYG"][...]); mxg = int(ds[0]["MXG"][...]); ny = int(ds[0]["ny"][...])
    t = float(ds[0]["t_array"][-1])/float(ds[0]["Omega_ci"][...])

    def get(name):
        parts = []
        for k, dd in enumerate(ds):
            v = dd[name]; a = np.asarray(v[-1, :, :, 0]) * (v.getncattr("conversion") if "conversion" in v.ncattrs() else 1.0)
            lo = 0 if k == 0 else myg; hi = a.shape[1] if k == len(ds) - 1 else a.shape[1] - myg
            parts.append(a[:, lo:hi])
        return np.concatenate(parts, axis=1)            # (nx, ny + 2 myg)
    N, Te, Ti, V = get("Nd+"), get("Te"), get("Td+"), get("Vd+")
    with netCDF4.Dataset(gridfile) as gf:
        R, Z = np.asarray(gf["Rxy"][:]), np.asarray(gf["Zxy"][:]); j1, j2 = int(gf["jyseps1_1"][...]), int(gf["jyseps2_2"][...]); ixs = int(gf["ixseps1"][...])
    xs = slice(mxg, R.shape[0] - mxg)
    inner = lambda a: a[xs, myg:-myg]                    # (nx - 2 mxg, ny)
    Rc, Zc, Nc, Tec, Tic = R[xs], Z[xs], inner(N), inner(Te), inner(Ti)
    psin_x = psin_RZ(Rc[:, -1], Zc[:, -1])               # psiN of each x index
    rows = np.arange(j1 + 1, j2 + 1)
    jm = rows[np.argmin(np.where(Rc[ixs - mxg, rows] > q["rmagx"], np.abs(Zc[ixs - mxg, rows] - Z_OMP), np.inf))]
    prof = {}
    xm = rsep_of_psin(psin_RZ(Rc[:, jm], Zc[:, jm])); cs = np.sqrt(qe*(Tec[:, jm] + Tic[:, jm])/(2*1.67262192e-27))
    prof["FHRP", "density"], prof["FHRP", "electron_temp"], prof["FHRP", "jsat"] = (xm, Nc[:, jm]), (xm, Tec[:, jm]), (xm, 0.5*qe*Nc[:, jm]*cs)
    sel = (np.abs(Rc - 0.9) < 0.012) & (Zc < -0.2) & (Zc > -0.70)
    xt = rsep_of_psin(psin_RZ(Rc[sel], Zc[sel]))
    prof["TS", "density"], prof["TS", "electron_temp"] = (xt, Nc[sel]), (xt, Tec[sel])
    xtar = rsep_of_psin(psin_x)
    for diag, jb in (("LFS-LP", ny + myg - 1), ("HFS-LP", myg - 1)):      # boundary between rows jb and jb + 1 of the guarded arrays
        b = lambda a: 0.5*(a[xs, jb] + a[xs, jb + 1])
        gam = b(N)*np.abs(b(V))
        prof[diag, "density"], prof[diag, "electron_temp"], prof[diag, "jsat"] = (xtar, b(N)), (xtar, b(Te)), (xtar, qe*gam)
        if diag == "LFS-LP":
            prof["LFS-IR", "q_parallel"] = (xtar, (GAMMA_E*b(Te) + GAMMA_I*b(Ti))*qe*gam + 13.6*qe*gam)
    isep = ixs - mxg
    print("%s (structured): t = %.1f ms | midplane row Z = %.3f m | either side of the separatrix: n_e %.2e / %.2e m^-3, Te %.1f / %.1f eV, Ti %.1f / %.1f eV" % (
        run, 1e3*t, Zc[isep, jm], Nc[isep - 1, jm], Nc[isep, jm], Tec[isep - 1, jm], Tec[isep, jm], Tic[isep - 1, jm], Tic[isep, jm]))
    for dd in ds: dd.close()
    return dict(struct=True, prof=prof, t=t)


def sim_profile(d, diag, obs):
    """(R^u - R^u_sep [cm], value) of the simulation for one observable, sorted"""
    if d.get("struct"):
        x, val = d["prof"][diag, obs]; o = np.argsort(x)
        return np.asarray(x)[o], np.asarray(val)[o]
    g, n = d["g"], d["n"]
    if diag == "FHRP":
        sel = np.flatnonzero((np.abs(g["Z"][:n] - Z_OMP) < 0.012) & (g["R"][:n] > q["rmagx"]))
        cs = np.sqrt(qe*(d["Te"] + d["Ti"])/(2*1.67262192e-27))
        val = {"density": d["Ne"], "electron_temp": d["Te"], "jsat": 0.5*qe*d["Ne"]*cs}[obs][sel]   # probe: J_sat = 0.5 e n c_s
        x = rsep_of_psin(d["psin"][sel])
    elif diag == "TS":
        sel = np.flatnonzero((np.abs(g["R"][:n] - 0.9) < 0.012) & (g["Z"][:n] < -0.2) & (g["Z"][:n] > -0.70))   # lowest Thomson point is at Z = -0.692
        val = {"density": d["Ne"], "electron_temp": d["Te"]}[obs][sel]
        x = rsep_of_psin(d["psin"][sel])
    else:
        fR, fZ = g["face_Rmid"], g["face_Zmid"]
        if diag in ("LFS-LP", "LFS-IR"):
            sel = np.flatnonzero(d["wall"] & (fZ < -0.749) & (fR > 0.69))       # floor, outboard of the private region's inner end
        else:
            sel = np.flatnonzero(d["wall"] & (fR < 0.6245) & (fZ < -0.2))       # inner wall below the midplane
        gh = d["R"][sel]
        psif = psin_RZ(fR[sel], fZ[sel])
        x = rsep_of_psin(psif)
        val = {"density": d["Ne"][gh], "electron_temp": d["Te"][gh], "jsat": qe*d["gam"][sel],
               "q_parallel": d["qpar"][sel] + 13.6*qe*d["gam"][sel]}[obs]
    o = np.argsort(x)
    return x[o], val[o]


def exp_profile(field, diag, obs):
    with netCDF4.Dataset(DATA + "1.experimental_data_TCV_%s_field.nc" % field) as d:
        o = d[diag]["observables"][obs]
        x, v, e = (np.asarray(o[k][:], float) for k in ("Rsep_omp", "value", "error"))
    ok = np.isfinite(x) & np.isfinite(v) & np.isfinite(e) & (e > 0)
    k = np.argsort(x[ok])
    return x[ok][k], v[ok][k], e[ok][k]


def distance(xs, vs, xe, ve, ee):
    inside = (xe > xs.min()) & (xe < xs.max())
    if inside.sum() < 2:
        return np.nan, 0
    s = np.interp(xe[inside], xs, vs)
    return float(np.sqrt(np.mean((s - ve[inside])**2/ee[inside]**2))), int(inside.sum())


if __name__ == "__main__":
    fig_path = sys.argv[1]
    runs = [a.split("=", 1) for a in sys.argv[2:]]
    data = {lab: (load_structured(*path.split("::")) if "::" in path else load(path)) for lab, path in runs}
    colors = ["#1f5fa8", "#2a8a57", "#7a4fb5", "#b8860b"]
    for lab, d in data.items():
        if d.get("struct"): continue
        dev = np.abs(d["psin"][:d["n"]] - d["grid_psin"][:d["n"]]).max()
        print("%s: t = %.1f ms; psiN of the grid file and of the geqdsk agree to %.1e; R_sep at Z=0: %.4f m" % (lab, 1e3*d["t"], dev, R_sep))

    panels = [("FHRP", "density", "midplane $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("FHRP", "electron_temp", "midplane $T_e$ [eV]", 1),
              ("TS", "density", "divertor entrance $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("TS", "electron_temp", "divertor entrance $T_e$ [eV]", 1),
              ("LFS-LP", "density", "outer target $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("LFS-LP", "electron_temp", "outer target $T_e$ [eV]", 1),
              ("LFS-LP", "jsat", "outer target $J_{sat}$ [kA m$^{-2}$]", 1e-3), ("LFS-IR", "q_parallel", "outer target $q_\\parallel$ [MW m$^{-2}$]", 1e-6),
              ("HFS-LP", "density", "inner target $n_e$ [10$^{18}$ m$^{-3}$]", 1e-18), ("HFS-LP", "electron_temp", "inner target $T_e$ [eV]", 1),
              ("HFS-LP", "jsat", "inner target $J_{sat}$ [kA m$^{-2}$]", 1e-3), ("FHRP", "jsat", "midplane $J_{sat}$ [kA m$^{-2}$]", 1e-3)]
    fig, axes = plt.subplots(3, 4, figsize=(19, 12.5))
    print("\n%-8s %-14s | %s | SOLPS-ITER standard (forward, reversed)" % ("diag", "observable", " | ".join("%s: d_j forward, reversed" % lab for lab in data)))
    summary = {lab: [] for lab in data}; solps_all = []
    for a, (diag, obs, title, scale) in zip(axes.ravel(), panels):
        xmin, xmax = np.inf, -np.inf
        exp = {}
        for field, c, mk in (("forward", "0.25", "o"), ("reversed", "0.6", "s")):
            xe, ve, ee = exp_profile(field, diag, obs); exp[field] = (xe, ve, ee)
            a.errorbar(xe, ve*scale, ee*scale, fmt=mk, ms=3, color=c, lw=0.8, capsize=1.5, label="experiment, %s field" % field)
            xmin, xmax = min(xmin, xe.min()), max(xmax, xe.max())
        row = []
        for (lab, d), c in zip(data.items(), colors):
            xs, vs = sim_profile(d, diag, obs)
            a.plot(xs, vs*scale, "-", color=c, lw=1.8, marker=".", ms=4, label="Hermes-3, " + lab)
            dj = [distance(xs, vs, *exp[f])[0] for f in ("forward", "reversed")]
            row.append("%6.2f %6.2f" % tuple(dj)); summary[lab].append(dj)
        solps_all.append(SOLPS[(diag, obs)])
        print("%-8s %-14s | %s | %.2f %.2f" % (diag, obs, " | ".join(row), *SOLPS[(diag, obs)]))
        pad = 0.15*(xmax - xmin); a.set_xlim(xmin - pad, xmax + pad); a.set_ylim(0, None)
        a.axvline(0, color="0.7", lw=0.8); a.set_title(title, fontsize=11)
    for a in axes[-1]: a.set_xlabel("$R^u - R^u_{sep}$ [cm]")
    axes[0, 0].legend(frameon=False, fontsize=8)
    for lab in data:
        m = np.nanmean(np.array(summary[lab]), axis=0)
        print("mean d_j over these %d observables: %s forward %.2f, reversed %.2f | SOLPS-ITER standard %.2f, %.2f" % (len(panels), lab, m[0], m[1], *np.mean(np.array(solps_all), axis=0)))
    fig.suptitle("TCV-X21: Hermes-3 against the experimental dataset (no drifts, so both field directions are shown)")
    fig.savefig(fig_path, dpi=105, bbox_inches="tight")
