"""Analytic test grids for the paper's numerical tests: an annulus of circular
flux surfaces (as tools/unstructured/annulus.py), made of skewed
quadrilaterals or of triangles.

    build(nr, ntheta, skew, jitter, triangles=False, psi_noise=0.0,
          outer="wall")

skew       cells of each ring are shifted along the flux surface by `skew`
           cell widths relative to the ring inside (a smooth, resolution-
           independent tilt of the faces that the field crosses)
jitter     random shift of every vertex along its flux surface (cell widths)
triangles  split every quadrilateral along a diagonal into two triangles
psi_noise  random error added to psi at every vertex, as a fraction of the
           psi range: the faces on flux surfaces then carry a small magnetic
           flux, as on a grid whose vertices are not snapped to the surfaces
outer      kind of the outer boundary: "wall" or "core"

Field: psi = PSI0 r^2 / 2 about (R0, 0), R Bt = RBT; B_pol is along +theta.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "unstructured"))
import flatgrid  # noqa: E402
from annulus import R0, R_IN, R_OUT, PSI0, RBT  # noqa: E402,F401


def build(nr=8, ntheta=32, skew=0.0, jitter=0.0, triangles=False, psi_noise=0.0,
          outer="wall", seed=1):
    rng = np.random.default_rng(seed)
    r = np.linspace(R_IN, R_OUT, nr + 1)
    dth = 2 * np.pi / ntheta
    th = (np.arange(ntheta)[None, :] + skew * np.arange(nr + 1)[:, None]
          + jitter * rng.uniform(-0.5, 0.5, (nr + 1, ntheta))) * dth
    vid = np.arange((nr + 1) * ntheta).reshape(nr + 1, ntheta)
    vx = (R0 + r[:, None] * np.cos(th)).ravel()
    vy = (r[:, None] * np.sin(th)).ravel()
    vpsi = (0.5 * PSI0 * r[:, None]**2 * np.ones_like(th)).ravel()
    if psi_noise > 0:
        span = 0.5 * PSI0 * (R_OUT**2 - R_IN**2)
        vpsi = vpsi + psi_noise * span * rng.uniform(-1, 1, vpsi.shape)
    vgrad = PSI0 * np.stack([vx - R0, vy], 1)
    jp = (np.arange(ntheta) + 1) % ntheta
    outer_kind = flatgrid.KIND_WALL if outer == "wall" else flatgrid.KIND_CORE

    cells, tubes, face_vx, face_cells, kind = [], [], [], [], []

    def face(v, c, k):
        face_vx.append(v)
        face_cells.append(c)
        kind.append(k)

    if not triangles:
        cid = np.arange(nr * ntheta).reshape(nr, ntheta)
        cells = [np.array([vid[i, j], vid[i + 1, j], vid[i + 1, jp[j]], vid[i, jp[j]]])
                 for i in range(nr) for j in range(ntheta)]
        tubes = [cid[i] for i in range(nr)]
        for i in range(nr):
            for j in range(ntheta):
                face((vid[i, jp[j]], vid[i + 1, jp[j]]), (cid[i, j], cid[i, jp[j]]), flatgrid.KIND_INTERIOR)
                if i + 1 < nr:
                    face((vid[i + 1, j], vid[i + 1, jp[j]]), (cid[i, j], cid[i + 1, j]), flatgrid.KIND_INTERIOR)
                else:
                    face((vid[i + 1, j], vid[i + 1, jp[j]]), (cid[i, j], -1), outer_kind)
        for j in range(ntheta):
            face((vid[0, j], vid[0, jp[j]]), (cid[0, j], -1), flatgrid.KIND_CORE)
    else:
        # quadrilateral (i, j) -> lower triangle A = (i,j),(i+1,j),(i+1,j+1)
        #                         upper triangle B = (i,j),(i+1,j+1),(i,j+1)
        A = lambda i, j: 2 * (i * ntheta + j)
        B = lambda i, j: 2 * (i * ntheta + j) + 1
        for i in range(nr):
            for j in range(ntheta):
                cells.append(np.array([vid[i, j], vid[i + 1, j], vid[i + 1, jp[j]]]))
                cells.append(np.array([vid[i, j], vid[i + 1, jp[j]], vid[i, jp[j]]]))
            tubes.append(np.array([c for j in range(ntheta) for c in (A(i, j), B(i, j))]))
        for i in range(nr):
            for j in range(ntheta):
                face((vid[i, j], vid[i + 1, jp[j]]), (A(i, j), B(i, j)), flatgrid.KIND_INTERIOR)         # diagonal
                face((vid[i, jp[j]], vid[i + 1, jp[j]]), (B(i, j), A(i, jp[j])), flatgrid.KIND_INTERIOR)  # radial
                if i + 1 < nr:   # flux surface between rings i and i+1
                    face((vid[i + 1, j], vid[i + 1, jp[j]]), (A(i, j), B(i + 1, j)), flatgrid.KIND_INTERIOR)
                else:
                    face((vid[i + 1, j], vid[i + 1, jp[j]]), (A(i, j), -1), outer_kind)
        for j in range(ntheta):
            face((vid[0, j], vid[0, jp[j]]), (B(0, j), -1), flatgrid.KIND_CORE)
    return flatgrid.build(vx, vy, vpsi, vgrad, RBT, cells, tubes, np.array(face_vx),
                          np.array(face_cells), np.array(kind))


def exact_unit_vectors(x):
    """Exact b (R, Z components of the unit vector along B, poloidal part) and
    the radial unit vector at points x (n, 2)"""
    rad = x - np.array([R0, 0.0])
    r = np.hypot(*rad.T)
    e_r = rad / r[:, None]
    e_th = np.stack([-e_r[:, 1], e_r[:, 0]], 1)
    Bp = PSI0 * r / x[:, 0]
    Bt = RBT / x[:, 0]
    B = np.hypot(Bp, Bt)
    return (Bp / B)[:, None] * e_th, e_r


if __name__ == "__main__":
    for tri in (False, True):
        g = build(8, 32, skew=0.4, triangles=tri)
        out = flatgrid.check(g, verbose=False)
        print("triangles" if tri else "quadrilaterals", {k: v for k, v in out.items()})
