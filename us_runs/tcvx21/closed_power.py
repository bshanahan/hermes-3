"""Heat crossing the closed flux surfaces of a TCV-X21 GOAT-grid run, by channel.

usage (from us_runs/): python tcvx21/closed_power.py <run dir>[@@<grid file>] [...]

The surfaces are the boundaries between neighbouring closed flux tubes of the
grid. Across each: heat conducted by electrons and ions (chi = 1 m^2/s),
carried by the anomalous particle diffusion (D = 0.2 m^2/s, 3/2 T per
particle), all with the operator the code uses (isotropic gradient, skew
correction, upwind coefficient), and carried by the E x B and magnetic drifts
(5/2 p v). Each drift flux is given as the code computes it (upwind) and
with the face value taken as the mean of the two cells (central); the
difference is the numerical diffusion of the upwinding.
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata, usops
import compare as c

qe = c.qe; CHI, D = 1.0, 0.2
rows = []
for arg in sys.argv[1:]:
    run, grid = (arg.split("@@") + ["grids/tcvx21_goat_drifts.grd.nc"])[:2]   # '<dir>@@<grid file>' for another grid
    r = usdata.Run(run, grid); g = r.g; n = r.n
    Ne = r("Nd+", True)[-1]; Te = r("Te", True)[-1]; Ti = r("Td+", True)[-1]
    phi = r("phi", True)[-1] if r.has("phi") else None
    tube = g["tube"]; L, R = g["face_L"], g["face_R"]
    closed = [t for t in np.unique(tube[:n]) if (g["prev"][:n][tube[:n] == t] >= 0).all() and (g["next"][:n][tube[:n] == t] >= 0).all()]
    psin = c.psin_RZ(g["R"][:n], g["Z"][:n])
    mean_psin = {t: psin[tube[:n] == t].mean() for t in closed}
    order = sorted(closed, key=lambda t: mean_psin[t])
    pe, pi = Ne*Te*qe, Ne*Ti*qe; P = pe + pi
    print("== %s (t = %.1f ms); %d closed tubes" % (run, 1e3*r.t[-1], len(order)))
    k = 2*np.pi/1e3
    heat = sum(1.5*r.volume_integral(r("P%s_src" % sp, True)[-1]) for sp in ("d+", "e") if r.has("P%s_src" % sp))/1e3
    line = "  heating sources %.1f kW" % heat
    for sp in ("d+", "e"):
        if r.has("P%s_core_heating_fraction" % sp):
            line += ", %s core heating fraction %.3f" % (sp, float(r.d["P%s_core_heating_fraction" % sp][-1]))
    fc = np.flatnonzero(g["face_kind"] == 1)
    if len(fc):
        Lf, Rf = L[fc], R[fc]
        def cflow(a, fld):                                    # plasma -> core ghost
            gr = usops.face_gradient(g, fld, isotropic=True, skew=True, limit=1.0, smooth=0.05, harmonic=False)[fc]
            return -gr*np.where(gr > 0, a[Rf], a[Lf])
        cond = (CHI*cflow(Ne, Te) + CHI*cflow(Ne, Ti))*qe
        conv = cflow(1.5*D*(Te + Ti), Ne)*qe
        part = cflow(D*np.ones_like(Ne), Ne)
        line += " | through the core boundary into the plasma: conduction %.1f kW, with particles %.1f kW, particles %.3e /s => total power in %.1f kW" % (
            -k*cond.sum(), -k*conv.sum(), -2*np.pi*part.sum(), heat - k*(cond.sum() + conv.sum()))
    print(line)
    print("  between psiN      cond e   cond i   anom. conv   ExB up/cen    mag up/cen    total(up)  [kW]  | phi spread on inner ring [V]")
    for a, b in zip(order[:-1], order[1:]):
        f = np.flatnonzero((g["face_kind"] == 0) & (((tube[L] == a) & (tube[R] == b)) | ((tube[L] == b) & (tube[R] == a))))
        out = np.where(tube[L[f]] == a, 1.0, -1.0)            # +1: L is the inner tube
        Lf, Rf = L[f], R[f]
        # the operator anomalous_diffusion uses: isotropic grad_perp, skew
        # correction, coefficient from the cell the flow leaves
        def flow(a, fld):          # L -> R
            gr = usops.face_gradient(g, fld, isotropic=True, skew=True, limit=1.0, smooth=0.05, harmonic=False)[f]
            return -gr*np.where(gr > 0, a[Rf], a[Lf])
        ce = (CHI*flow(Ne, Te)*qe*out).sum(); ci = (CHI*flow(Ne, Ti)*qe*out).sum()
        conv = (flow(1.5*D*(Te + Ti), Ne)*qe*out).sum()      # 3/2 T per particle, as the code
        row = {"run": run, "psin": (mean_psin[a], mean_psin[b]), "cond": k*(ce + ci), "conv": k*conv}
        line = "  %.3f-%.3f   %7.1f  %7.1f    %7.1f" % (mean_psin[a], mean_psin[b], k*ce, k*ci, k*conv)
        tot = ce + ci + conv
        if phi is not None:
            fE = usops.vertex_difference(g, phi)[f]*g["face_drift"][f]          # L -> R
            up = lambda flow, p: np.where(flow > 0, p[Lf], p[Rf])
            cen = lambda p: 0.5*(p[Lf] + p[Rf])
            Tfi = 0.5*(Ti[Lf] + Ti[Rf]); Tfe = 0.5*(Te[Lf] + Te[Rf])
            fDi = Tfi*g["face_curl"][f]; fDe = -Tfe*g["face_curl"][f]
            exb_u = (2.5*up(fE, P)*fE*out).sum(); exb_c = (2.5*cen(P)*fE*out).sum()
            mag_u = (2.5*(up(fDi, pi)*fDi + up(fDe, pe)*fDe)*out).sum(); mag_c = (2.5*(cen(pi)*fDi + cen(pe)*fDe)*out).sum()
            tot += exb_u + mag_u
            ring = np.flatnonzero(tube[:n] == a)
            line += "   %6.1f/%6.1f  %6.1f/%6.1f" % (k*exb_u, k*exb_c, k*mag_u, k*mag_c)
            row.update(exb_u=k*exb_u, exb_c=k*exb_c, mag_u=k*mag_u, mag_c=k*mag_c)
            spread = "%.0f..%.0f" % (phi[ring].min(), phi[ring].max())
        else:
            line += "   %13s  %13s" % ("-", "-"); spread = "-"
        row["total"] = k*tot
        print(line + "   %7.1f        | %s" % (k*tot, spread))
        rows.append(row)
