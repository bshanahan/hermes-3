"""Analytic test grid: an annulus of circular flux surfaces, optionally skewed.

usage: python annulus.py out.grd.nc [nr ntheta skew jitter]

The flux surfaces are circles about (R0, 0); psi = psi0 * r^2 / 2. The inner
circle is the core boundary and the outer one the wall, so the field does not
cross any boundary. With `skew` the cells of each ring are shifted along the
flux surface by `skew` cell widths per ring relative to the ring inside, and
`jitter` moves every vertex along its flux surface by a random fraction of a
cell width. Both leave the flux surfaces exact and make the line between
neighbouring cell centres oblique to their common face.
"""
import sys
import numpy as np
import flatgrid

R0, R_IN, R_OUT, PSI0, RBT = 1.5, 0.2, 0.4, 0.5, 3.0


def build(nr=8, ntheta=32, skew=0.0, jitter=0.0, seed=1):
    rng = np.random.default_rng(seed)
    r = np.linspace(R_IN, R_OUT, nr + 1)
    dth = 2 * np.pi / ntheta
    th = (np.arange(ntheta)[None, :] + skew * np.arange(nr + 1)[:, None]
          + jitter * rng.uniform(-0.5, 0.5, (nr + 1, ntheta))) * dth
    vid = np.arange((nr + 1) * ntheta).reshape(nr + 1, ntheta)
    vx = (R0 + r[:, None] * np.cos(th)).ravel()
    vy = (r[:, None] * np.sin(th)).ravel()
    vpsi = (0.5 * PSI0 * r[:, None]**2 * np.ones_like(th)).ravel()
    vgrad = PSI0 * np.stack([vx - R0, vy], 1)

    cid = np.arange(nr * ntheta).reshape(nr, ntheta)
    jp = (np.arange(ntheta) + 1) % ntheta
    cells = [np.array([vid[i, j], vid[i + 1, j], vid[i + 1, jp[j]], vid[i, jp[j]]])
             for i in range(nr) for j in range(ntheta)]
    tubes = [cid[i] for i in range(nr)]
    face_vx, face_cells, kind = [], [], []
    for i in range(nr):
        for j in range(ntheta):
            # face between (i, j) and (i, j+1), crossed by the field
            face_vx.append((vid[i, jp[j]], vid[i + 1, jp[j]]))
            face_cells.append((cid[i, j], cid[i, jp[j]]))
            kind.append(flatgrid.KIND_INTERIOR)
            # face on the flux surface outside (i, j)
            face_vx.append((vid[i + 1, j], vid[i + 1, jp[j]]))
            if i + 1 < nr:
                face_cells.append((cid[i, j], cid[i + 1, j]))
                kind.append(flatgrid.KIND_INTERIOR)
            else:
                face_cells.append((cid[i, j], -1))
                kind.append(flatgrid.KIND_WALL)
    for j in range(ntheta):
        face_vx.append((vid[0, j], vid[0, jp[j]]))
        face_cells.append((cid[0, j], -1))
        kind.append(flatgrid.KIND_CORE)
    return flatgrid.build(vx, vy, vpsi, vgrad, RBT, cells, tubes, np.array(face_vx),
                          np.array(face_cells), np.array(kind))


if __name__ == "__main__":
    a = sys.argv[2:]
    g = build(*(int(v) for v in a[:2]), *(float(v) for v in a[2:4]))
    flatgrid.check(g)
    flatgrid.write(g, sys.argv[1], blocks=12)
    print("wrote", sys.argv[1])
