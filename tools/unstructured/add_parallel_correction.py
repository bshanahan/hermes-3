"""Add the coefficients of the corrected parallel conduction (us_face_kparn,
us_face_vpar; see flatgrid.build) to a flat-cell grid file written before they
existed. Cells, their order and all other variables stay as they are, so restart
files remain valid.

usage: python add_parallel_correction.py <grid file> [<grid file> ...]

The poloidal field direction at a face is rebuilt from what the file holds: its
component along the face normal from the flux through the face, the size of the
other component from |B_pol|, and its sign from the two cells beside the face.
"""
import sys

import netCDF4
import numpy as np

KIND_INTERIOR = 0


def coefficients(d):
    xc = np.stack([d["us_R"][:], d["us_Z"][:]], 1)                # storage order
    P = np.stack([d["us_vertex_R"][:], d["us_vertex_Z"][:]], 1)
    L, R = np.asarray(d["us_face_L"][:]), np.asarray(d["us_face_R"][:])
    v1, v2 = np.asarray(d["us_face_v1"][:]), np.asarray(d["us_face_v2"][:])
    kind = np.asarray(d["us_face_kind"][:])
    flux, Bf, Bp = (np.asarray(d["us_face_" + k][:]) for k in ("flux", "B", "Bp"))
    farea, dn = np.asarray(d["us_face_area"][:]), np.asarray(d["us_face_dn"][:])
    nrm = np.stack([d["us_face_nR"][:], d["us_face_nZ"][:]], 1)
    interior = kind == KIND_INTERIOR
    t = P[v2] - P[v1]
    tau = t / np.hypot(*t.T)[:, None]
    dvec = xc[R] - xc[L]
    det = dvec[:, 0] * t[:, 1] - dvec[:, 1] * t[:, 0]
    beta = np.stack([-dvec[:, 1], dvec[:, 0]], 1) / det[:, None]
    bn = flux / (farea * Bf)
    bt_size = np.sqrt(np.maximum((Bp / Bf)**2 - bn**2, 0.0))
    bcell = np.stack([d["us_bR"][:], d["us_bZ"][:]], 1)
    Rr = np.where(interior, R, L)                                 # ghost cells carry no direction
    sign = np.sign(np.einsum("ij,ij->i", bcell[L] + bcell[Rr], tau))
    b = bn[:, None] * nrm + (sign * bt_size)[:, None] * tau       # in-plane part of the unit vector b
    # the two-point form stays between two quadrilaterals (see flatgrid.build)
    nreal = int(d["us_ncells"][...])
    l_of_s = -np.ones(len(xc), int)
    s_of_l = np.asarray(d["us_logical"][:])
    l_of_s[s_of_l] = np.arange(len(s_of_l))
    nvert = np.full(len(xc), 4)
    nvert[s_of_l[:nreal]] = (np.asarray(d["us_cell_vertices"][:]) >= 0).sum(1)
    skewed = interior & ((nvert[L] != 4) | (nvert[R] != 4))
    kpar = np.asarray(d["us_face_kpar"][:])
    kparn = np.where(skewed, (flux / Bf) * bn / dn, kpar)
    vpar = np.where(skewed, (flux / Bf) * np.einsum("ij,ij->i", b, beta), 0.0)
    return kparn, vpar


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with netCDF4.Dataset(path, "a") as d:
            if "us_face_sptr" not in d.variables:
                print(path, ": no vertex stencils, skipped")
                continue
            kparn, vpar = coefficients(d)
            for name, val in (("us_face_kparn", kparn), ("us_face_vpar", vpar)):
                v = d[name] if name in d.variables else d.createVariable(name, "f8", ("us_face",))
                v[:] = val
            kpar = np.asarray(d["us_face_kpar"][:])
            print("%s: %d faces carry parallel conduction, %d of them (beside a cell that is not a "
                  "quadrilateral) now use the diamond gradient" % (path, (kpar > 0).sum(), (vpar != 0).sum()))
