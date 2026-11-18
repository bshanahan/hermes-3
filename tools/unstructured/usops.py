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


def vertex_difference(g, f):
    """Value at vertex 2 minus value at vertex 1 of every face, from the cell stencils"""
    ptr = g["face_sptr"]
    face = np.repeat(np.arange(len(ptr) - 1), np.diff(ptr))
    out = np.zeros(len(ptr) - 1)
    np.add.at(out, face, g["face_sw"] * f[g["face_scell"]])
    return out


def face_gradient(g, f, isotropic=False, skew=False, limit=1.0, smooth=0.0):
    """Transmissibility times the difference of f across every face.

    skew: add the vertex correction, kept below `limit` times the two-point
    term with the same smooth bound as the C++ (limit=None: no bound).
    smooth: this fraction of the mean size of f counts as a difference between
    the two cells when the bound is formed (unstructured:skew_smooth).
    """
    L, R = g["face_L"], g["face_R"]
    trans = g["face_tiso"] if isotropic else g["face_tperp"]
    two_point = trans * (f[R] - f[L])
    if not skew:
        return two_point
    corr = (g["face_viso"] if isotropic else g["face_vperp"]) * vertex_difference(g, f)
    if limit is None:
        return two_point + corr
    floor_diff = smooth * 0.5 * (np.abs(f[L]) + np.abs(f[R])) * trans
    bound = limit * np.sqrt(two_point**2 + floor_diff**2)
    return two_point + corr * bound / np.sqrt(bound**2 + corr**2 + 1e-300)


def div_a_grad_perp_upwind(g, a, f, isotropic=False, skew=False, limit=1.0, smooth=0.0):
    L, R = g["face_L"], g["face_R"]
    grad = face_gradient(g, f, isotropic, skew, limit, smooth)
    flow_into_L = grad * np.where(grad > 0, a[R], a[L])
    return _scatter(g, flow_into_L)


def wall_loss(g, a, f, inv_length):
    """Loss through wall faces with a decay length; returns the (negative) source"""
    wall = g["face_kind"] == 2
    L = g["face_L"]
    bn_pol = np.minimum(1.0, np.abs(g["face_flux"]) / (g["face_area"] * np.maximum(g["face_Bp"], 1e-30)))
    out = np.where(wall, a[L] * f[L] * g["face_area"] * np.sqrt(1 - bn_pol**2) * inv_length, 0.0)
    res = np.zeros(len(g["vol"]))
    np.add.at(res, L, -out)
    res /= g["vol"]
    res[g["nreal"]:] = 0.0
    return res


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


def div_a_grad_perp(g, a, f, skew=False, limit=1.0, smooth=0.0):
    """Div(a Grad_perp f) with the mean coefficient at each face, full grad_perp"""
    L, R = g["face_L"], g["face_R"]
    flow_into_L = 0.5 * (a[L] + a[R]) * face_gradient(g, f, True, skew, limit, smooth)
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
