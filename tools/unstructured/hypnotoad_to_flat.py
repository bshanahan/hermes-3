"""Recast a structured hypnotoad grid as a flat-cell grid with the same cells.

This is for verification: the flat-cell (unstructured) operators in Hermes-3
should then reproduce a structured run on the original grid.

Single-null and lower/upper single-null topologies with x guard cells in the
file are supported.

usage: python hypnotoad_to_flat.py bout.grd.nc equilibrium.geqdsk out.grd.nc [--metric]

With --metric the cell volumes, face areas and transmissibilities are taken
from the structured metric (J, dx, dy, g11, g22) instead of the cell polygons.
"""
import sys
import numpy as np
import netCDF4
from scipy.interpolate import RectBivariateSpline
import flatgrid


def convert(gridfile, geqdsk, mxg=2, metric=False):
    from hypnotoad.geqdsk._geqdsk import read as read_geqdsk

    d = netCDF4.Dataset(gridfile)
    nx, ny = len(d.dimensions["x"]), len(d.dimensions["y"])
    ix1 = int(d["ixseps1"][...])
    ix2 = int(d["ixseps2"][...])
    j11, j21, j12, j22 = (int(d[k][...]) for k in ("jyseps1_1", "jyseps2_1", "jyseps1_2", "jyseps2_2"))
    if ix2 < nx - mxg or j21 != j12:
        raise ValueError("only single-null grids are supported")

    corner = {
        "ll": (d["Rxy_corners"][:], d["Zxy_corners"][:]),
        "lr": (d["Rxy_lower_right_corners"][:], d["Zxy_lower_right_corners"][:]),
        "ul": (d["Rxy_upper_left_corners"][:], d["Zxy_upper_left_corners"][:]),
        "ur": (d["Rxy_upper_right_corners"][:], d["Zxy_upper_right_corners"][:]),
    }
    psi_xlow = np.asarray(d["psixy_xlow"][:])
    psixy = np.asarray(d["psixy"][:])
    # psi on the x-high face of the last cell: extrapolate
    psi_xhigh = np.empty_like(psi_xlow)
    psi_xhigh[:-1] = psi_xlow[1:]
    psi_xhigh[-1] = 2 * psixy[-1] - psi_xlow[-1]
    rbt = float(np.mean(np.asarray(d["Btxy"][:]) * np.asarray(d["Rxy"][:])))

    # Unique vertices, matched by position
    vid = {}
    vR, vZ, vpsi = [], [], []

    def vertex(name, i, j):
        R, Z = float(corner[name][0][i, j]), float(corner[name][1][i, j])
        key = (round(R, 9), round(Z, 9))
        if key not in vid:
            vid[key] = len(vR)
            vR.append(R)
            vZ.append(Z)
            vpsi.append(psi_xlow[i, j] if name in ("ll", "ul") else psi_xhigh[i, j])
        return vid[key]

    xs, xe = mxg, nx - mxg - 1
    cell_id = -np.ones((nx, ny), int)
    cells = []
    for i in range(xs, xe + 1):
        for j in range(ny):
            cell_id[i, j] = len(cells)
            cells.append(np.array([vertex("ll", i, j), vertex("lr", i, j), vertex("ur", i, j), vertex("ul", i, j)]))

    # Faces, matched by their vertex pair
    fid = {}
    face_vx, face_cells, face_kind = [], [], []

    def face(a, b, cell, kind):
        key = (min(a, b), max(a, b))
        if key in fid:
            f = fid[key]
            assert face_cells[f][1] < 0
            face_cells[f][1] = cell
            face_kind[f] = flatgrid.KIND_INTERIOR
        else:
            fid[key] = len(face_vx)
            face_vx.append([a, b])
            face_cells.append([cell, -1])
            face_kind.append(kind)

    def in_core(j):
        return j11 < j <= j22

    for i in range(xs, xe + 1):
        for j in range(ny):
            c = cell_id[i, j]
            ll, lr, ur, ul = cells[c]
            # Kind used if the face turns out to have only this cell
            face(ll, ul, c, flatgrid.KIND_CORE if (i == xs and in_core(j)) else flatgrid.KIND_WALL)  # x low
            face(lr, ur, c, flatgrid.KIND_WALL)  # x high
            face(ll, lr, c, flatgrid.KIND_WALL)  # y low
            face(ul, ur, c, flatgrid.KIND_WALL)  # y high

    # Flux tubes
    tubes = []
    for i in range(xs, xe + 1):
        if i < ix1:
            tubes.append(cell_id[i, j11 + 1:j22 + 1])  # closed
            tubes.append(np.concatenate([cell_id[i, :j11 + 1], cell_id[i, j22 + 1:]]))  # private flux
        else:
            tubes.append(cell_id[i, :])

    vR, vZ, vpsi = np.array(vR), np.array(vZ), np.array(vpsi)

    # Gradient of psi at the vertices from the equilibrium
    with open(geqdsk) as f:
        q = read_geqdsk(f)
    r = q["rleft"] + q["rdim"] * np.arange(q["nx"]) / (q["nx"] - 1)
    z = q["zmid"] - q["zdim"] / 2 + q["zdim"] * np.arange(q["ny"]) / (q["ny"] - 1)
    spl = RectBivariateSpline(r, z, np.asarray(q["psi"]), kx=3, ky=3)
    vgrad = np.stack([spl(vR, vZ, dx=1, grid=False), spl(vR, vZ, dy=1, grid=False)], 1)
    # The grid's psi may differ from the geqdsk's by sign and offset; match the sign
    sample = spl(vR, vZ, grid=False)
    slope = np.polyfit(sample, vpsi, 1)[0]
    vgrad *= np.sign(slope)
    print("grid psi = %.6f * geqdsk psi + const (expected +-1)" % slope)

    g = flatgrid.build(vR, vZ, vpsi, vgrad, rbt, cells, tubes, np.array(face_vx), np.array(face_cells),
                       np.array(face_kind))
    g["structured_index"] = np.argwhere(cell_id >= 0)  # (i, j) of each cell in the original order
    order = np.concatenate([np.asarray(t) for t in tubes])
    g["structured_index"] = g["structured_index"][order]
    if metric:
        use_structured_metric(g, d)
    return g


def use_structured_metric(g, d):
    """Replace volumes, areas and transmissibilities by the values the
    structured finite-volume operators of BOUT++/Hermes-3 use on this grid
    (J, dx, dy, g11, g22, g_22), so that the flat-cell operators can be
    compared with the structured ones without geometric differences."""
    n = g["nreal"]
    xi, yj = g["structured_index"].T
    J, dx, dy, g11, g22, g_22 = (np.asarray(d[k][:]) for k in ("J", "dx", "dy", "g11", "g22", "g_22"))
    g["vol"][:n] = (J * dx * dy)[xi, yj]
    L, R = g["face_L"], g["face_R"]
    il, jl = xi[L], yj[L]
    real = R < n
    Rr = np.where(real, R, L)  # boundary faces: use the cell itself for the far side
    ir, jr = xi[Rr], yj[Rr]
    hy_l = np.sqrt(g_22[il, jl])
    hy_r = np.sqrt(g_22[ir, jr])
    dxf = 0.5 * (dx[il, jl] + dx[ir, jr])
    # Faces crossed by the field are the constant-y faces. (Radial faces have
    # a round-off level flux, so test the angle rather than for zero.)
    tube = np.abs(g["face_flux"] / (g["face_area"] * g["face_Bp"])) > 0.5
    radial = ~tube
    # Area for parallel flows, J dx / sqrt(g_22) at the face
    apar = (J[il, jl] + J[ir, jr]) / (hy_l + hy_r) * dxf
    g["face_flux"] = np.where(tube, np.sign(g["face_flux"]) * apar * g["face_B"], 0.0)
    # Parallel distance between cell centres
    dpar = 0.5 * (dy[il, jl] * hy_l + dy[ir, jr] * hy_r)
    g["face_kpar"] = np.where(tube & real, apar / dpar, 0.0)
    # Radial faces: (J g11) / dx * dy
    trad = (J[il, jl] * g11[il, jl] + J[ir, jr] * g11[ir, jr]) / (dx[il, jl] + dx[ir, jr]) * dy[il, jl]
    # Poloidal part of grad_perp through constant-y faces: J (g22 - 1/g_22) / dy * dx
    perp_y = 0.5 * (J[il, jl] * (g22[il, jl] - 1 / g_22[il, jl]) + J[ir, jr] * (g22[ir, jr] - 1 / g_22[ir, jr]))
    tpol = perp_y / (0.5 * (dy[il, jl] + dy[ir, jr])) * dxf
    wall = g["face_kind"] == flatgrid.KIND_WALL
    core = g["face_kind"] == flatgrid.KIND_CORE
    # A core ghost is at the face: half the centre-to-centre distance
    g["face_tperp"] = np.where(radial & ~wall, trad * np.where(core, 2.0, 1.0), 0.0)
    g["face_tiso"] = np.where(wall, 0.0, np.where(radial, g["face_tperp"], tpol))


if __name__ == "__main__":
    grid = convert(sys.argv[1], sys.argv[2], metric="--metric" in sys.argv)
    print("%d cells + %d ghost cells, %d faces" % (grid["nreal"], grid["nghost"], grid["nfaces"]))
    flatgrid.check(grid, strict="--metric" not in sys.argv)
    ij = grid["structured_index"]
    pad = np.zeros(grid["nreal"] + grid["nghost"])
    xi, yj = pad.copy(), pad.copy()
    xi[:grid["nreal"]], yj[:grid["nreal"]] = ij[:, 0], ij[:, 1]
    flatgrid.write(grid, sys.argv[3], extra={"structured_x": xi, "structured_y": yj})
    print("wrote", sys.argv[3])
