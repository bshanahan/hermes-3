"""Size of the cross-field flux error on skewed cells, with and without the vertex correction.

usage: python skew_report.py traduit.out.b2us

For a field that varies linearly in (R, Z) the flow through every face is
known exactly. The two-point flow uses the two cell centres only; the
corrected flow adds the vertex term. Errors are given relative to
D * |grad f| * (face area), for a gradient along the radial direction, along
the poloidal field and (isotropic operator) along and across the face normal.
"""
import sys
import numpy as np
import goat_to_bout


def errors(g, which, grad):
    """relative error of the two-point and corrected flows for unit gradients `grad` (nf, 2)"""
    nreal = g["nreal"]
    xc = np.stack([g["R"], g["Z"]], 1)
    L, R = g["face_L"], g["face_R"]
    d = xc[R] - xc[L]
    t = g["face_t" + which] * np.einsum("ij,ij->i", grad, d)
    # vertex term: sum over the stencil of sw * f(cell), f = grad . x
    corr = np.zeros(len(L))
    ptr, cell, sw = g["face_sptr"], g["face_scell"], g["face_sw"]
    for f in range(len(L)):
        k = slice(ptr[f], ptr[f + 1])
        corr[f] = g["face_v" + which][f] * np.dot(sw[k], xc[cell[k]] @ grad[f])
    exact = np.einsum("ij,ij->i", g["nvec_" + which], grad) * g["face_area"]
    return np.abs(t - exact) / g["face_area"], np.abs(t + corr - exact) / g["face_area"]


def report(g):
    nrm = np.stack([g["face_nR"], g["face_nZ"]], 1)
    P = np.stack([g["vertex_R"], g["vertex_Z"]], 1)
    bn = g["face_flux"] / (g["face_area"] * g["face_B"])
    bn_pol = np.clip(g["face_flux"] / (g["face_area"] * np.maximum(g["face_Bp"], 1e-30)), -1, 1)
    # unit vector along the poloidal field at each face, from the face data: flux = area * Bp * (bhat . n)
    tng = P[g["face_v2"]] - P[g["face_v1"]]; tng /= np.hypot(*tng.T)[:, None]
    sgn = np.sign(np.einsum("ij,ij->i", g["bhat_face"], tng)); sgn[sgn == 0] = 1
    bhat = g["bhat_face"]
    g["nvec_perp"] = nrm - bn_pol[:, None] * bhat
    g["nvec_iso"] = nrm - bn[:, None] * bhat * (g["face_Bp"] / g["face_B"])[:, None]
    er = np.stack([bhat[:, 1], -bhat[:, 0]], 1)
    interior = g["face_kind"] == 0
    tube = g["face_side_L"] != 0
    out = {}
    for which, sel, cases in (("perp", interior & ~tube, (("radial gradient", er), ("poloidal gradient", bhat))),
                              ("perp", interior & tube, (("radial gradient", er),)),
                              ("iso", interior, (("gradient along the face normal", nrm), ("gradient along the face", tng)))):
        for label, grad in cases:
            e2, ec = errors(g, which, grad)
            name = "%s operator, %s faces, %s" % ({"perp": "radial", "iso": "isotropic"}[which],
                                                   "all interior" if which == "iso" else ("between-tube" if not tube[sel][0] else "along-tube"), label)
            q = lambda e: tuple(np.round(np.percentile(e[sel], [50, 90, 99, 100]), 3))
            out[name] = (int(sel.sum()), q(e2), q(ec))
    vk = g["vertex_kind"]
    print("%d faces (%d interior), %d vertices (%d interior, %d core, %d wall, %d without a linear fit)" % (len(tube), interior.sum(), len(vk), (vk % 10 == 0).sum(), (vk % 10 == 1).sum(), (vk % 10 == 2).sum(), (vk >= 10).sum()))
    print("relative flow error for a linear field: median, 90th, 99th percentile, maximum")
    for k, (n, a, b) in out.items():
        print("  %-72s n=%5d  two-point %s   corrected %s" % (k, n, a, b))
    return out




if __name__ == "__main__":
    import flatgrid
    from b2us import read_b2us as read_b2us
    raw = read_b2us(sys.argv[1])
    g = goat_to_bout.convert(sys.argv[1])
    # poloidal field direction at faces, as flatgrid.build computes it
    vgrad = np.stack([raw["vdpsidx"], raw["vdpsidy"]], 1) / (2 * np.pi)
    gp = 0.5 * (vgrad[g["face_v1"]] + vgrad[g["face_v2"]])
    bpol = np.stack([-gp[:, 1], gp[:, 0]], 1)
    g["bhat_face"] = bpol / np.hypot(*bpol.T)[:, None]
    report(g)
