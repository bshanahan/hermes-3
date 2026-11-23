"""Does parallel conduction leak across flux surfaces? (cf. Guenter et al. 2005)

usage: python pollution.py out.json

Steady anisotropic conduction  div(k_par b b.grad T + k_perp grad_perp T) + S = 0
on the skewed annulus, with the flat-cell operators (two-point parallel flow,
corrected cross-field flow without the bound, so that the problem is linear).
The manufactured solution T(r) is constant on flux surfaces, so b.grad T = 0
and the exact solution does not depend on k_par / k_perp. Any growth of the
error with k_par / k_perp is numerical transport across the flux surfaces.

Cases: vertices exactly on the flux surfaces ("snapped") or with a psi error
of 3e-4 of the psi range (as GOAT places them before snapping), on
quadrilaterals and on triangles.
"""
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annuli  # noqa: E402
import flatgrid  # noqa: E402

RATIOS = [0.0] + [10.0**k for k in range(0, 11)]
W = annuli.R_OUT - annuli.R_IN


def T_exact(x):
    r = np.hypot(x[:, 0] - annuli.R0, x[:, 1])
    return 2.0 + np.cos(np.pi * (r - annuli.R_IN) / W)


def grad_T(x):
    rad = x - np.array([annuli.R0, 0.0])
    r = np.hypot(*rad.T)
    dT = -np.pi / W * np.sin(np.pi * (r - annuli.R_IN) / W)
    return (dT / r)[:, None] * rad


def cell_sources(g, kperp):
    """Exact integral of S over each cell (per radian): minus the outward
    exact flow k_perp grad T . n through its faces, by Gauss quadrature"""
    xg, wg = np.polynomial.legendre.leggauss(8)
    P = np.stack([g["vertex_R"], g["vertex_Z"]], 1)
    v1, v2 = g["face_v1"], g["face_v2"]
    flow = np.zeros(len(v1))          # exact flow from L to R through each face, along n (L->R)
    n = np.stack([g["face_nR"], g["face_nZ"]], 1)
    for s, w in zip(0.5 * (xg + 1), 0.5 * wg):
        x = P[v1] + s * (P[v2] - P[v1])
        flow += w * (-kperp * np.einsum("ij,ij->i", grad_T(x), n)) * x[:, 0]
    flow *= g["face_length"]
    src = np.zeros(len(g["vol"]))
    np.add.at(src, g["face_L"], flow)          # what leaves L must be made up by the source
    np.add.at(src, g["face_R"], -flow)
    return src[:g["nreal"]]


def beta_vectors(g):
    """beta_f: in the poloidal plane, beta.d = 0 and beta.t = 1 (t from vertex 1 to vertex 2)"""
    P = np.stack([g["vertex_R"], g["vertex_Z"]], 1)
    t = P[g["face_v2"]] - P[g["face_v1"]]
    d = np.stack([g["R"][g["face_R"]] - g["R"][g["face_L"]], g["Z"][g["face_R"]] - g["Z"][g["face_L"]]], 1)
    det = d[:, 0] * t[:, 1] - d[:, 1] * t[:, 0]
    return np.stack([-d[:, 1], d[:, 0]], 1) / det[:, None]


def solve(g, kpar, kperp=1.0, parallel="two_point"):
    nreal = g["nreal"]
    ntot = len(g["vol"])
    L, R, kind = g["face_L"], g["face_R"], g["face_kind"]
    x = np.stack([g["R"], g["Z"]], 1)
    known = np.zeros(ntot)
    core = kind == flatgrid.KIND_CORE
    known[R[core]] = T_exact(x[R[core]])        # ghost cells sit at the face mid-points
    rows, cols, vals = [], [], []
    rhs = -cell_sources(g, kperp)

    def add(i, j, v):                            # coefficient of T_j in the inflow to cell i
        if i >= nreal:
            return
        if j < nreal:
            rows.append(i); cols.append(j); vals.append(v)
        else:
            rhs[i] -= v * known[j]

    interior = kind == flatgrid.KIND_INTERIOR
    c = kperp * g["face_tiso"] + np.where(interior, kpar * g["face_kpar"], 0.0)
    vpar = np.zeros(len(L))
    if parallel == "diamond":
        # parallel flow from the diamond gradient: (Phi/B) b.G, G as for the cross-field correction
        xm = np.stack([g["face_Rmid"], g["face_Zmid"]], 1)
        b, _ = annuli.exact_unit_vectors(xm)
        nrm = np.stack([g["face_nR"], g["face_nZ"]], 1)
        fb = g["face_flux"] / g["face_B"]
        bn = np.einsum("ij,ij->i", b, nrm)
        c = kperp * g["face_tiso"] + np.where(interior, kpar * fb * bn / g["face_dn"], 0.0)
        vpar = np.where(interior, kpar * fb * np.einsum("ij,ij->i", b, beta_vectors(g)), 0.0)
    for f in range(len(L)):
        if c[f] == 0 and g["face_viso"][f] == 0 and vpar[f] == 0:
            continue
        i, j = L[f], R[f]
        # inflow to i = c (T_j - T_i) + kperp viso (T_v2 - T_v1)
        add(i, j, c[f]); add(i, i, -c[f])
        if interior[f]:
            add(j, i, c[f]); add(j, j, -c[f])
        if g["face_viso"][f] != 0 or vpar[f] != 0:
            for k in range(g["face_sptr"][f], g["face_sptr"][f + 1]):
                v = (kperp * g["face_viso"][f] + vpar[f]) * g["face_sw"][k]
                add(i, g["face_scell"][k], v)
                if interior[f]:
                    add(j, g["face_scell"][k], -v)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(nreal, nreal))
    T = spla.spsolve(A.tocsc(), rhs)
    ex = T_exact(x[:nreal])
    w = g["vol"][:nreal]
    return float(np.sqrt(np.sum(w * (T - ex)**2) / np.sum(w)) / (ex.max() - ex.min()))


if __name__ == "__main__":
    out = {"ratios": RATIOS, "cases": {}}
    for tri in (False, True):
        for noise in (0.0, 3e-4):
            for n in (16, 32):
                g = annuli.build(n, 4 * n, skew=0.4, triangles=tri, psi_noise=noise)
                key = "%s, %s, n=%d" % ("triangles" if tri else "quadrilaterals",
                                        "psi error 3e-4" if noise else "snapped", n)
                for par in ("two_point", "diamond"):
                    if par == "diamond" and noise:
                        continue
                    errs = [solve(g, k, parallel=par) for k in RATIOS]
                    out["cases"][key + ("" if par == "two_point" else ", diamond parallel")] = errs
                    print(key, par, " ".join("%.1e" % e for e in errs), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
