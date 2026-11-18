"""Convert a GOAT grid (traduit.out.b2us) into a flat-cell grid file for Hermes-3.

usage: python goat_to_bout.py traduit.out.b2us out.grd.nc [--tube-order]
"""
import sys
import numpy as np
from b2us import read_b2us
import flatgrid


def snap_flux_surfaces(vpsi, face_vx, aligned):
    """Give all vertices joined by flux-surface faces exactly the same psi.

    GOAT places these vertices on a flux surface only to a tolerance (a few
    1e-4 of the local gradient times the face length). The magnetic flux
    through a face is the psi difference of its vertices, so the residue
    lets parallel transport leak through faces that should carry none. At
    keV temperatures that leak is far larger than the physical cross-field
    transport. Vertices connected by aligned faces form one flux surface
    (the separatrix branches join at the X-point); each such group is set to
    its mean psi. Fluxes of every cell still sum to zero exactly.
    """
    parent = np.arange(len(vpsi))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for (a, b), flag in zip(face_vx, aligned):
        if flag:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    roots = np.array([find(i) for i in range(len(vpsi))])
    out = vpsi.copy()
    for root in np.unique(roots):
        members = roots == root
        if members.sum() > 1:
            out[members] = vpsi[members].mean()
    return out, float(np.abs(out - vpsi).max())


def convert(b2us_path):
    g = read_b2us(b2us_path)
    twopi = 2 * np.pi
    # GOAT stores psi in Wb and its gradient in Wb/m, and 2 pi R Bt
    vpsi = g["vpsi"] / twopi
    vgrad = np.stack([g["vdpsidx"], g["vdpsidy"]], 1) / twopi
    rbt = g["vffbz"][0] / twopi

    vpsi, spread = snap_flux_surfaces(vpsi, g["fcVx"], g["fcAligned"])
    print("vertices moved onto their flux surface: largest psi change %.2e Wb/rad "
          "(psi range %.2f)" % (spread, np.ptp(vpsi)))

    face_cells = -np.ones((g["nFc"], 2), int)
    for c, faces in enumerate(g["cell_faces"]):
        for f in faces:
            face_cells[f, 0 if face_cells[f, 0] < 0 else 1] = c
    assert np.all(face_cells[:, 0] >= 0)
    # GOAT labels: 0 interior, < 0 core boundary, > 0 vessel structure
    kind = np.where(g["fcLbl"] == 0, flatgrid.KIND_INTERIOR,
                    np.where(g["fcLbl"] < 0, flatgrid.KIND_CORE, flatgrid.KIND_WALL))
    assert np.all((kind == flatgrid.KIND_INTERIOR) == (face_cells[:, 1] >= 0))
    return flatgrid.build(g["vx"], g["vy"], vpsi, vgrad, rbt, g["cells"], g["tubes"], g["fcVx"],
                          face_cells, kind, g["fcLbl"])


if __name__ == "__main__":
    grid = convert(sys.argv[1])
    print("%d cells + %d ghost cells, %d faces" % (grid["nreal"], grid["nghost"], grid["nfaces"]))
    flatgrid.check(grid)
    flatgrid.write(grid, sys.argv[2], order="tube" if "--tube-order" in sys.argv else "angle")
    print("wrote", sys.argv[2])
