"""A GOAT grid restricted to a range of poloidal flux: the flux tubes of
traduit.out.b2us whose mean psi lies inside the range, as a flat-cell grid file.
The faces towards the removed tubes become wall faces. Closed flux tubes are
always kept.

This gives an unstructured grid of the region that an orthogonal (structured)
grid covers, to compare the two kinds of grid on the same domain.

usage: python goat_cut.py traduit.out.b2us out.grd.nc <psi min> <psi max>   [Wb/rad, as in the grid file]
"""
import sys
import numpy as np
from b2us import read_b2us
import flatgrid
from goat_to_bout import snap_flux_surfaces


def convert(b2us_path, psi_min, psi_max):
    g = read_b2us(b2us_path)
    twopi = 2 * np.pi
    vpsi = g["vpsi"] / twopi
    vgrad = np.stack([g["vdpsidx"], g["vdpsidy"]], 1) / twopi
    rbt = g["vffbz"][0] / twopi
    vpsi, _ = snap_flux_surfaces(vpsi, g["fcVx"], g["fcAligned"])

    ncell = len(g["cells"])
    face_cells = -np.ones((g["nFc"], 2), int)
    for c, faces in enumerate(g["cell_faces"]):
        for f in faces:
            face_cells[f, 0 if face_cells[f, 0] < 0 else 1] = c
    kind = np.where(g["fcLbl"] == 0, flatgrid.KIND_INTERIOR,
                    np.where(g["fcLbl"] < 0, flatgrid.KIND_CORE, flatgrid.KIND_WALL))

    # Tubes to keep: closed ones (those that touch no wall face and whose cells all have
    # two tube neighbours are the ones inside the separatrix: psi below it) and the
    # open ones inside the range
    cpsi = np.array([vpsi[np.asarray(c, int)].mean() for c in g["cells"]])
    wall_cell = np.zeros(ncell, bool)
    wall_cell[face_cells[kind == flatgrid.KIND_WALL, 0]] = True
    core_cell = np.zeros(ncell, bool)
    core_cell[face_cells[kind == flatgrid.KIND_CORE, 0]] = True
    psi_core = cpsi[core_cell].mean()
    tubes, removed = [], 0
    for t in g["tubes"]:
        t = np.asarray(t, int)
        p = cpsi[t].mean()
        closed = not wall_cell[t].any()
        if closed or psi_min < p < psi_max:
            tubes.append(t)
        else:
            removed += 1
    print("%d of %d flux tubes kept (%d removed); core boundary cells at psi %.5f" % (
        len(tubes), len(g["tubes"]), removed, psi_core))
    keep = np.zeros(ncell, bool)
    keep[np.concatenate(tubes)] = True
    new_cell = -np.ones(ncell, int)
    new_cell[keep] = np.arange(keep.sum())

    # Faces with at least one kept cell; a face that loses a neighbour is a wall
    fc, fk, fv, fl = [], [], [], []
    cut = 0
    for f in range(g["nFc"]):
        a, b = face_cells[f]
        ka, kb = keep[a], (b >= 0 and keep[b])
        if not (ka or kb):
            continue
        if ka and kb:
            fc.append([new_cell[a], new_cell[b]]); fk.append(kind[f])
        elif b < 0:
            fc.append([new_cell[a], -1]); fk.append(kind[f])
        else:
            fc.append([new_cell[a if ka else b], -1]); fk.append(flatgrid.KIND_WALL)
            cut += 1
        fv.append(g["fcVx"][f]); fl.append(g["fcLbl"][f] if (b < 0 or (ka and kb)) else 1)
    print("%d faces towards removed tubes are now wall faces" % cut)
    fv = np.asarray(fv, int)

    # Vertices in use
    cells = [np.asarray(g["cells"][c], int) for c in np.flatnonzero(keep)]
    used = np.zeros(len(vpsi), bool)
    for c in cells:
        used[c] = True
    new_v = -np.ones(len(vpsi), int)
    new_v[used] = np.arange(used.sum())
    cells = [new_v[c] for c in cells]
    fv = new_v[fv]
    assert fv.min() >= 0
    tubes = [new_cell[t] for t in tubes]
    return flatgrid.build(g["vx"][used], g["vy"][used], vpsi[used], vgrad[used], rbt, cells, tubes, fv,
                          np.asarray(fc, int), np.asarray(fk), np.asarray(fl))


if __name__ == "__main__":
    grid = convert(sys.argv[1], float(sys.argv[3]), float(sys.argv[4]))
    print("%d cells + %d ghost cells, %d faces" % (grid["nreal"], grid["nghost"], grid["nfaces"]))
    flatgrid.check(grid)
    flatgrid.write(grid, sys.argv[2], order="angle")
    print("wrote", sys.argv[2])
