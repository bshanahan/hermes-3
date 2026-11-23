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


def face_gradient(g, f, isotropic=False, skew=False, limit=1.0, smooth=1.0, harmonic=True,
                  power=2.0, scale="positive"):
    """Transmissibility times the difference of f across every face.

    skew: add the vertex correction C, bounded as in the C++
    (Unstructured::correctedGradient): C / (1 + |C/M|^power)^(1/power) with
        M = limit * sqrt(A^2 + (smooth * trans * h)^2),
    A the two-point term and h the harmonic (or arithmetic) mean of |f| in the
    two cells (limit=None: no bound). The defaults are those of the code
    (unstructured:skew_limit, skew_smooth, skew_harmonic_mean, skew_bound_power);
    smooth=0.05, harmonic=False is the bound of the first version.
    scale: what kind of field f is (Unstructured::Scale): "positive" (harmonic
    mean), "log" (h = 1) or "signed" (root mean square); only with harmonic=True.
    """
    L, R = g["face_L"], g["face_R"]
    trans = g["face_tiso"] if isotropic else g["face_tperp"]
    two_point = trans * (f[R] - f[L])
    if not skew:
        return two_point
    corr = (g["face_viso"] if isotropic else g["face_vperp"]) * vertex_difference(g, f)
    if limit is None:
        return two_point + corr
    a, b = np.abs(f[L]), np.abs(f[R])
    if not harmonic:
        h = 0.5 * (a + b)
    elif scale == "positive":
        h = 2 * a * b / (a + b + 1e-300)
    elif scale == "log":
        h = np.ones_like(a)
    elif scale == "signed":
        h = np.sqrt(0.5 * (f[L]**2 + f[R]**2))
    else:
        raise ValueError(scale)
    bound = limit * np.sqrt(two_point**2 + (smooth * h * trans)**2)
    if power == 2.0:
        return two_point + corr * bound / np.sqrt(bound**2 + corr**2 + 1e-300)
    r = np.abs(corr) / (bound + 1e-300)
    return two_point + corr / (1 + r**power)**(1.0 / power)


def div_a_grad_perp_upwind(g, a, f, isotropic=False, skew=False, limit=1.0, smooth=1.0,
                           harmonic=True, power=2.0):
    L, R = g["face_L"], g["face_R"]
    grad = face_gradient(g, f, isotropic, skew, limit, smooth, harmonic, power)
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


def div_a_grad_perp(g, a, f, skew=False, limit=1.0, smooth=1.0, harmonic=True, power=2.0):
    """Div(a Grad_perp f) with the mean coefficient at each face, full grad_perp"""
    L, R = g["face_L"], g["face_R"]
    flow_into_L = 0.5 * (a[L] + a[R]) * face_gradient(g, f, True, skew, limit, smooth, harmonic,
                                                      power)
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


def flow_exb(g, phi, through_wall=True, phi_width=0.0):
    """E x B volume flow through every face, from L to R. Through wall faces only
    the flow into the wall is kept (none if `through_wall` is false); with
    `phi_width` [V] that cut-off is smooth, as in Hermes-3 (1e-3 * Tnorm)."""
    flow = g["face_drift"] * vertex_difference(g, phi)
    wall = g["face_kind"] == 2
    width = phi_width * np.abs(g["face_drift"][wall])
    flow[wall] = 0.5 * (flow[wall] + np.sqrt(flow[wall]**2 + width**2)) if through_wall else 0.0
    return flow


def cell_gradient(g, f):
    """Green-Gauss gradient (R, Z) in every real cell, as Unstructured::cellGradient"""
    L, R, kind = g["face_L"], g["face_R"], g["face_kind"]
    n = len(g["vol"])
    bnd = kind != 0
    w = np.where(bnd, f[R], 0.5 * (f[L] + f[R])) * g["face_length"]
    sR, sZ = np.zeros(n), np.zeros(n)
    np.add.at(sR, L, w * g["face_nR"])
    np.add.at(sZ, L, w * g["face_nZ"])
    np.add.at(sR, R[~bnd], -(w * g["face_nR"])[~bnd])
    np.add.at(sZ, R[~bnd], -(w * g["face_nZ"])[~bnd])
    real = np.arange(n) < g["nreal"]
    gR, gZ = np.zeros(n), np.zeros(n)
    gR[real] = sR[real] / g["area"][real]
    gZ[real] = sZ[real] / g["area"][real]
    return gR, gZ


def limited_face_value(fc, fd, two_g_dot_d):
    """Face value from the upwind cell, TVD with a smoothed van Albada limiter
    (as limitedFaceValue in src/unstructured.cxx)"""
    b = fd - fc
    a = two_g_dot_d - b
    scale = 1e-6 * (np.abs(fc) + np.abs(fd)) + 1e-300
    ab = a * b
    positive = 0.5 * (ab + np.sqrt(ab**2 + scale**4))
    return fc + 0.5 * positive * (a + b) / (a**2 + b**2 + scale**2 + 1e-300)


def drift_face_values(g, f, flow, phi_width=0.0, second_order=True, n=None, n_blend=0.0):
    """Value of f carried by a drift flow through every face, as Unstructured::Div_f_flow.
    With a density n and n_blend > 0 the second-order value is blended with the
    upwind one by n_f/(n_f + n_blend), n_f the harmonic mean of the two cells
    (Unstructured::secondOrderWeight)."""
    L, R, kind = g["face_L"], g["face_R"], g["face_kind"]
    core = kind == 1
    fcore = (g["face_area"][core] * f[R[core]]).sum() / g["face_area"][core].sum() if core.any() else 0.0
    width = phi_width * np.abs(g["face_drift"]) if "face_drift" in g else np.zeros_like(flow)
    wL = 0.5 * (1 + flow / np.sqrt(flow**2 + width**2 + 1e-300))
    if second_order:
        gR, gZ = cell_gradient(g, f)
        dR, dZ = g["R"][R] - g["R"][L], g["Z"][R] - g["Z"][L]
        fromL = limited_face_value(f[L], f[R], 2 * (gR[L] * dR + gZ[L] * dZ))
        fromR = limited_face_value(f[R], f[L], -2 * (gR[R] * dR + gZ[R] * dZ))
        inner = wL * fromL + (1 - wL) * fromR
        if n is not None and n_blend > 0:
            nL, nR = np.maximum(n[L], 0), np.maximum(n[R], 0)
            nf = 2 * nL * nR / (nL + nR + 1e-300)
            w = nf / (nf + n_blend)
            inner = (wL * f[L] + (1 - wL) * f[R]) + w * (inner - (wL * f[L] + (1 - wL) * f[R]))
    else:
        inner = wL * f[L] + (1 - wL) * f[R]
    return np.where(core, fcore, np.where(kind == 2, f[L], inner))


def div_f_flow(g, f, flow, phi_width=0.0, second_order=True):
    """Div(f v) for a volume flow through each face (Unstructured::Div_f_flow):
    limited second-order (or upwind) values on interior faces, the
    area-weighted mean of the boundary values on core faces, the value of the
    cell on wall faces."""
    L, R = g["face_L"], g["face_R"]
    value = drift_face_values(g, f, flow, phi_width, second_order)
    out = np.zeros(len(g["vol"]))
    np.add.at(out, L, flow * value)
    np.add.at(out, R, -flow * value)
    out /= g["vol"]
    out[g["nreal"]:] = 0.0
    return out
