"""Accuracy of the flat-cell face flows on skewed analytic annuli (numpy
reference operators, tools/unstructured/usops.py, which the C++ reproduces).

usage: python op_convergence.py out.json

For skewed quadrilaterals and for triangles at several resolutions, compares
the flow through every interior face with the exact flow of an analytic field:
  cross-field   S n.(I - bb).grad u      two-point, corrected, corrected+bound
  coefficient   a S n.(I - bb).grad u    a from the upwind cell (plasma) or
                                         the mean of the two cells (neutrals)
  parallel      (Phi/B) b.grad u         two-point along the field
Errors are RMS over interior faces whose vertices are both interior, as a
fraction of the largest |grad u| (or |b.grad u|) on those faces.
Also: the largest error for linear fields (exactness).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annuli  # noqa: E402
import flatgrid  # noqa: E402
import usops  # noqa: E402

SKEW = 0.4
RES = (4, 8, 16, 32, 64)
BOUND = dict(limit=1.0, smooth=1.0, harmonic=True, power=2.0)      # the code's defaults (v18)
BOUND_V1 = dict(limit=1.0, smooth=0.05, harmonic=False, power=2.0)  # bound of the first version


def fields(name, x):
    """value and gradient (R, Z) of a test field at points x"""
    R, Z = x[:, 0], x[:, 1]
    if name == "smooth":
        k = 9.0
        u = 3.0 + np.sin(k * R) * np.cos(k * Z)
        g = np.stack([k * np.cos(k * R) * np.cos(k * Z), -k * np.sin(k * R) * np.sin(k * Z)], 1)
        return u, g
    if name == "steep":       # edge-like: e-folding length 3 cm across the 20 cm annulus
        lam = 0.03
        rad = x - np.array([annuli.R0, 0.0])
        r = np.hypot(*rad.T)
        th = np.arctan2(rad[:, 1], rad[:, 0])
        e = np.exp(-(r - annuli.R_IN) / lam)
        a = 1 + 0.2 * np.cos(th)
        u = e * a
        du_dr = -u / lam
        du_dth = e * (-0.2 * np.sin(th))
        e_r = rad / r[:, None]
        e_th = np.stack([-e_r[:, 1], e_r[:, 0]], 1)
        return u, du_dr[:, None] * e_r + (du_dth / r)[:, None] * e_th
    raise ValueError(name)


def face_data(g):
    L, R = g["face_L"], g["face_R"]
    xm = np.stack([g["face_Rmid"], g["face_Zmid"]], 1)
    n = np.stack([g["face_nR"], g["face_nZ"]], 1)
    b, _ = annuli.exact_unit_vectors(xm)
    sel = ((g["face_kind"] == flatgrid.KIND_INTERIOR)
           & (g["vertex_kind"][g["face_v1"]] == 0) & (g["vertex_kind"][g["face_v2"]] == 0))
    return L, R, xm, n, b, sel


def run(tri, n):
    t0 = time.time()
    g = annuli.build(n, 4 * n, skew=SKEW, triangles=tri)
    L, R, xm, nrm, b, sel = face_data(g)
    x = np.stack([g["R"], g["Z"]], 1)
    h = np.sqrt(np.mean(g["area"][:g["nreal"]]))
    out = {"n": n, "cells": int(g["nreal"]), "h": float(h)}
    for name in ("smooth", "steep"):
        u, _ = fields(name, x)
        _, gm = fields(name, xm)
        # exact flows: S n.(I - bb).grad u   (b has only its poloidal part in the plane)
        ngrad = np.einsum("ij,ij->i", nrm, gm)
        bn = np.einsum("ij,ij->i", nrm, b)
        bg = np.einsum("ij,ij->i", b, gm)
        exact = g["face_area"] * (ngrad - bn * bg)
        scale = np.abs(np.hypot(*gm.T))[sel].max()
        err = lambda a: float(np.sqrt(np.mean(((a - exact) / g["face_area"])[sel]**2)) / scale)
        two = usops.face_gradient(g, u, True)
        cor = usops.face_gradient(g, u, True, skew=True, limit=None)
        lim = usops.face_gradient(g, u, True, skew=True, **BOUND)
        lim1 = usops.face_gradient(g, u, True, skew=True, **BOUND_V1)
        out[name] = {"two_point": err(two), "corrected": err(cor), "bounded": err(lim),
                     "bounded_v1": err(lim1)}
        # fraction of faces where the bound changes the correction by more than 10%
        corr = cor - two
        out[name]["bound_active"] = float(np.mean((np.abs(lim - two) < 0.9 * np.abs(corr))[sel & (np.abs(corr) > 1e-14 * np.abs(exact).max())]))
        # coefficient a = u at the face: upwind (plasma) or mean (neutrals), corrected gradient
        a_exact = fields(name, xm)[0]
        up = np.where(cor > 0, u[R], u[L]) * cor
        mean = 0.5 * (u[L] + u[R]) * cor
        ex_a = a_exact * exact
        scale_a = np.abs(a_exact * np.hypot(*gm.T))[sel].max()
        err_a = lambda a: float(np.sqrt(np.mean(((a - ex_a) / g["face_area"])[sel]**2)) / scale_a)
        out[name]["coef_upwind"] = err_a(up)
        out[name]["coef_mean"] = err_a(mean)
        # parallel: (Phi/B) b.grad u, on faces the field crosses
        par = g["face_kpar"] * (u[R] - u[L])
        fb = g["face_flux"] / g["face_B"]
        crossed = sel & (np.abs(fb) > 0)
        ex_p = fb * bg
        scale_p = np.abs(bg)[crossed].max()
        out[name]["parallel"] = float(np.sqrt(np.mean(((par - ex_p) / np.abs(fb))[crossed]**2)) / scale_p)
    # exactness for linear fields u = 3 + G.x (offset so that the bound is formed as for a positive field)
    worst = {"two_point": 0.0, "corrected": 0.0, "bounded": 0.0, "bounded_v1": 0.0}
    for G in (np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.6, -0.8])):
        u = 3.0 + x @ G
        exact = g["face_area"] * (nrm @ G - np.einsum("ij,ij->i", nrm, b) * (b @ G))
        for key, val in (("two_point", usops.face_gradient(g, u, True)),
                         ("corrected", usops.face_gradient(g, u, True, skew=True, limit=None)),
                         ("bounded", usops.face_gradient(g, u, True, skew=True, **BOUND)),
                         ("bounded_v1", usops.face_gradient(g, u, True, skew=True, **BOUND_V1))):
            worst[key] = max(worst[key], float((np.abs(val - exact) / g["face_area"])[sel].max()))
    out["linear_max_error"] = worst
    out["seconds"] = time.time() - t0
    return out


if __name__ == "__main__":
    results = {"quadrilaterals": [], "triangles": []}
    for tri in (False, True):
        for n in RES:
            r = run(tri, n)
            results["triangles" if tri else "quadrilaterals"].append(r)
            print("triangles" if tri else "quads", n, r["cells"], "%.1fs" % r["seconds"],
                  {k: (v if not isinstance(v, dict) else {kk: "%.2e" % vv for kk, vv in v.items()})
                   for k, v in r.items() if k in ("smooth", "steep", "linear_max_error")}, flush=True)
    json.dump(results, open(sys.argv[1], "w"), indent=1)
