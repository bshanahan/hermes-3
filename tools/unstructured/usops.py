"""numpy reference implementation of the flat-cell operators (SI units).

Independent of the C++ in src/unstructured.cxx; used to check it.
Fields have one value per cell, ghosts included.
"""
import numpy as np


def _scatter(g, flow, boundary_into_L_only=True):
    """Divergence from face flows L->R: result = -(net outflow)/vol"""
    n = len(g["vol"])
    out = np.zeros(n)
    np.add.at(out, g["face_L"], flow)
    np.add.at(out, g["face_R"], -flow)
    out /= g["vol"]
    out[g["nreal"]:] = 0.0
    return out


def div_a_grad_perp_upwind(g, a, f, isotropic=False):
    L, R = g["face_L"], g["face_R"]
    grad = (g["face_tiso"] if isotropic else g["face_tperp"]) * (f[R] - f[L])
    flow_into_L = grad * np.where(grad > 0, a[R], a[L])
    return _scatter(g, flow_into_L)


def div_par_k_grad_par(g, K, f):
    L, R = g["face_L"], g["face_R"]
    flow_into_L = 0.5 * (K[L] + K[R]) * g["face_kpar"] * (f[R] - f[L])
    flow_into_L[g["face_kind"] != 0] = 0.0
    return _scatter(g, flow_into_L)


def div_par(g, f):
    L, R = g["face_L"], g["face_R"]
    fb = g["face_flux"] / g["face_B"]
    face = np.where(g["face_kind"] != 0, f[R], 0.5 * (f[L] + f[R]))
    return -_scatter(g, -face * fb) * -1 if False else _scatter(g, face * fb)


def grad_par(g, f):
    L, R = g["face_L"], g["face_R"]
    fb = g["face_flux"] / g["face_B"]
    n = len(g["vol"])
    out = np.zeros(n)
    half = np.where(g["face_kind"] != 0, 1.0, 0.5) * (f[R] - f[L]) * fb
    np.add.at(out, L, half)
    interior = g["face_kind"] == 0
    np.add.at(out, R[interior], half[interior])
    out /= g["vol"]
    out[g["nreal"]:] = 0.0
    return out


def div_a_grad_perp(g, a, f):
    """Div(a Grad_perp f) with the mean coefficient at each face, full grad_perp"""
    L, R = g["face_L"], g["face_R"]
    flow_into_L = 0.5 * (a[L] + a[R]) * g["face_tiso"] * (f[R] - f[L])
    return _scatter(g, flow_into_L)


def grad_perp_sq(g, f):
    """|Grad_perp f|^2 at real cells from a Green-Gauss gradient in the poloidal plane"""
    L, R = g["face_L"], g["face_R"]
    n = len(g["vol"])
    boundary = g["face_kind"] != 0
    face = np.where(boundary, f[R], 0.5 * (f[L] + f[R])) * g["face_length"]
    gR, gZ = np.zeros(n), np.zeros(n)
    np.add.at(gR, L, face * g["face_nR"])
    np.add.at(gZ, L, face * g["face_nZ"])
    np.add.at(gR, R[~boundary], -(face * g["face_nR"])[~boundary])
    np.add.at(gZ, R[~boundary], -(face * g["face_nZ"])[~boundary])
    gR, gZ = gR / g["area"], gZ / g["area"]
    par = g["bR"] * gR + g["bZ"] * gZ
    out = np.maximum(gR**2 + gZ**2 - par**2, 0.0)
    out[g["nreal"]:] = 0.0
    return out
