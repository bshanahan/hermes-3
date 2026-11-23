"""Flat-cell ("unstructured") grid files for Hermes-3.

The grid is a list of cells stored along BOUT++'s y index (nx = 1, MXG = 0).
Real cells come first, ordered flux tube by flux tube, followed by one ghost
cell per boundary face. A ghost cell holds the value AT its boundary face.

All geometry is in SI units and per radian of toroidal angle:

  us_vol          cell volume / 2 pi                         [m^3]
  us_face_area    face area / 2 pi                           [m^2]
  us_face_flux    poloidal magnetic flux through the face,
                  from cell L to cell R                      [Wb / rad]
  us_face_B       |B| at the face                            [T]
  us_face_kpar    parallel transmissibility:
                  conductive flow L->R = -kappa * kpar * (T_R - T_L)   [m]
  us_face_tperp   radial transmissibility:
                  diffusive flow L->R = -D * tperp * (n_R - n_L)       [m]
  us_face_tiso    as tperp, but for all of grad_perp (not only radial) [m]
  us_face_vperp, us_face_viso
                  skew correction: the flow gains -D * vperp * (n_v2 - n_v1),
                  where n_v is the value at a vertex of the face          [m]
  us_face_sptr, us_face_scell, us_face_sw
                  (n_v2 - n_v1) as a weighted sum over cells, one CSR row per face
  us_face_drift   E x B: the volume flow L->R through the face is
                  drift * (phi_v2 - phi_v1) = (R Bt / B^2) * (phi_v2 - phi_v1),
                  signed for the direction of the face normal. On a wall face
                  L->R is out of the plasma                      [m / T]
  us_face_curl    flux of curl(b / B) through the face from L to R, per radian:
                  the difference of R Bt / B^2 between the two vertices. The
                  magnetic (grad-B and curvature) drift carries
                  (p / q) * curl through the face. Zero on wall faces, and the
                  values of every cell without a wall face sum to zero  [m / T]
  us_area         poloidal area of the cell                  [m^2]
  us_face_length  poloidal length of the face                [m]
  us_bR, us_bZ    poloidal components of the unit vector along B at the cell

The flux through a face is the difference of psi at its two vertices, so the
fluxes of every cell sum to zero exactly (div B = 0 holds discretely) and
faces lying on a flux surface carry no parallel flow.
"""
import numpy as np
import netCDF4

KIND_INTERIOR, KIND_CORE, KIND_WALL, KIND_PADDING = 0, 1, 2, 3


def build(vx, vy, vpsi, vgradpsi, rbt, cells, tubes, face_vx, face_cells, face_kind,
          face_label=None):
    """Assemble the flat grid.

    vx, vy          vertex coordinates [m]
    vpsi            poloidal flux at vertices [Wb/rad]
    vgradpsi        (dpsi/dR, dpsi/dZ) at vertices [Wb/rad/m], shape (nv, 2)
    rbt             R * B_toroidal [T m]
    cells           list of vertex-index arrays, ordered around each polygon
    tubes           list of cell-index arrays, each ordered along the field
    face_vx         (nf, 2) vertex indices
    face_cells      (nf, 2) cell indices, -1 where there is no second cell
    face_kind       (nf,) KIND_* ; boundary faces must have face_cells[:, 1] == -1
    """
    P = np.stack([vx, vy], 1)
    ncell_old = len(cells)

    # New cell order: tube by tube
    order = np.concatenate([np.asarray(t, int) for t in tubes])
    assert len(order) == ncell_old and len(set(order.tolist())) == ncell_old
    new_of_old = np.empty(ncell_old, int)
    new_of_old[order] = np.arange(ncell_old)
    cells = [np.asarray(cells[o], int) for o in order]
    nreal = ncell_old

    # Cell centroid, area, volume per radian (integral of R dA)
    cR, cZ, area, vol = (np.zeros(nreal) for _ in range(4))
    for i, c in enumerate(cells):
        x, y = vx[c], vy[c]
        xn, yn = np.roll(x, -1), np.roll(y, -1)
        cr = x * yn - xn * y
        a = 0.5 * cr.sum()
        cR[i] = ((x + xn) * cr).sum() / (6 * a)
        cZ[i] = ((y + yn) * cr).sum() / (6 * a)
        area[i] = abs(a)
        vol[i] = abs(((x * x + x * xn + xn * xn) * (yn - y)).sum() / 6.0)

    # Tube neighbours
    prev = -np.ones(nreal, int)
    nxt = -np.ones(nreal, int)
    tube_id = np.zeros(nreal, int)
    pairs = set()
    fc = np.asarray(face_cells, int).copy()
    real = fc >= 0
    fc[real] = new_of_old[fc[real]]
    for a_, b_ in fc[fc[:, 1] >= 0]:
        pairs.add((min(a_, b_), max(a_, b_)))
    start = 0
    nclosed = 0
    for it, t in enumerate(tubes):
        n = len(t)
        idx = np.arange(start, start + n)
        tube_id[idx] = it
        prev[idx[1:]] = idx[:-1]
        nxt[idx[:-1]] = idx[1:]
        if n > 2 and (idx[0], idx[-1]) in pairs:  # closed flux tube
            prev[idx[0]], nxt[idx[-1]] = idx[-1], idx[0]
            nclosed += 1
        start += n

    # Faces
    nf = len(face_vx)
    face_vx = np.asarray(face_vx, int)
    face_kind = np.asarray(face_kind, int)
    v1, v2 = face_vx.T
    mid = 0.5 * (P[v1] + P[v2])
    t = P[v2] - P[v1]
    length = np.hypot(*t.T)
    nrm = np.stack([t[:, 1], -t[:, 0]], 1) / length[:, None]
    flux = -(vpsi[v2] - vpsi[v1])  # through the face along nrm
    farea = length * mid[:, 0]
    gpsi = 0.5 * (vgradpsi[v1] + vgradpsi[v2])
    bpol = np.stack([-gpsi[:, 1], gpsi[:, 0]], 1) / mid[:, [0]]  # (B_R, B_Z)
    Bp = np.hypot(*bpol.T)
    Bt = rbt / mid[:, 0]
    Bf = np.hypot(Bp, Bt)

    L = fc[:, 0].copy()
    R = fc[:, 1].copy()
    bnd = np.flatnonzero(R < 0)
    assert np.all(face_kind[bnd] != KIND_INTERIOR) and np.all(face_kind[R >= 0] == KIND_INTERIOR)
    nghost = len(bnd)
    R[bnd] = nreal + np.arange(nghost)
    ntot = nreal + nghost

    # Positions of all cells; a ghost sits at its face midpoint
    xc = np.zeros((ntot, 2))
    xc[:nreal] = np.stack([cR, cZ], 1)
    xc[R[bnd]] = mid[bnd]
    d = xc[R] - xc[L]
    flip = np.einsum("ij,ij->i", nrm, d) < 0
    nrm[flip] *= -1
    flux[flip] *= -1
    dn = np.einsum("ij,ij->i", nrm, d)

    # Drifts. With (R, phi, Z) right-handed and the normal (t_Z, -t_R) of the
    # face from vertex 1 to vertex 2, the poloidal E x B velocity
    # (Bt / B^2) e_phi x grad(phi) carries the volume R Bt / B^2 * (phi_2 - phi_1)
    # through the face, and Stokes' theorem gives the flux of curl(b / B) as
    # -(R Bt / B^2 at vertex 2 - at vertex 1).
    orient = np.where(flip, -1.0, 1.0)
    vB2 = ((vgradpsi**2).sum(1) + rbt**2) / vx**2          # B^2 at the vertices
    fdrift = orient * rbt / Bf**2
    fcurl = -orient * rbt * (1.0 / vB2[v2] - 1.0 / vB2[v1])
    fcurl[face_kind == KIND_WALL] = 0.0

    bn = flux / (farea * Bf)  # b . n
    bn_pol = np.clip(flux / (farea * np.maximum(Bp, 1e-30)), -1, 1)  # b_pol-hat . n
    kpar = np.maximum((flux / Bf) * np.einsum("ij,ij->i", bpol, d) / (Bf * np.einsum("ij,ij->i", d, d)), 0.0)
    tperp = farea * (1 - bn_pol**2) / dn
    tiso = farea * (1 - bn**2) / dn
    tperp[face_kind == KIND_WALL] = 0.0
    tiso[face_kind == KIND_WALL] = 0.0

    # Skew correction. The two-point flow above takes the gradient from the two
    # cell centres only, which is exact when the line joining them is normal to
    # the face. In general the gradient G in the diamond (L, v1, R, v2) follows
    # from   G.d = f_R - f_L,   G.t = f_v2 - f_v1   (t = x_v2 - x_v1):
    #   G = (f_R - f_L) n / dn + (f_v2 - f_v1) beta,   beta = (-d_y, d_x) / (d x t)
    # so the centre term is unchanged and the vertex term is the correction.
    det = d[:, 0] * t[:, 1] - d[:, 1] * t[:, 0]
    beta = np.stack([-d[:, 1], d[:, 0]], 1) / det[:, None]
    n_r = nrm - bn_pol[:, None] * bpol / np.maximum(Bp, 1e-30)[:, None]   # radial part of n
    n_iso = nrm - bn[:, None] * bpol / Bf[:, None]                        # n without its parallel part
    vperp = farea * np.einsum("ij,ij->i", n_r, beta)
    viso = farea * np.einsum("ij,ij->i", n_iso, beta)
    vperp[face_kind == KIND_WALL] = 0.0
    viso[face_kind == KIND_WALL] = 0.0
    sptr, scell, sw, vertex_kind = _vertex_stencils(P, xc, cells, v1, v2, L, R, face_kind, nreal)

    side_L = np.where(R == nxt[np.minimum(L, nreal - 1)], 1, np.where(R == prev[np.minimum(L, nreal - 1)], -1, 0))
    Rr = np.minimum(R, nreal - 1)
    side_R = np.where(R >= nreal, 0, np.where(L == nxt[Rr], 1, np.where(L == prev[Rr], -1, 0)))

    # Upstream cell for boundary extrapolation: the tube neighbour of L that
    # lies behind L as seen from the face, if the field crosses the face.
    up = -np.ones(nf, int)
    for f in bnd:
        i = L[f]
        best, bestval = -1, 0.3 * np.hypot(*d[f])
        for j in (prev[i], nxt[i]):
            if j >= 0:
                val = np.dot(xc[i] - xc[j], nrm[f])
                if val > bestval:
                    best, bestval = j, val
        if abs(bn_pol[f]) > 0.05:
            up[f] = best

    # Cell quantities including ghosts
    def cellfield(real_vals, ghost_vals):
        a = np.zeros(ntot)
        a[:nreal] = real_vals
        a[R[bnd]] = ghost_vals
        return a

    # B at cell centres from the vertices of each cell
    cgp = np.array([vgradpsi[c].mean(0) for c in cells])
    cBp = np.hypot(*cgp.T) / cR
    cBt = rbt / cR
    cpsi = np.array([vpsi[c].mean() for c in cells])

    cbvec = np.stack([-cgp[:, 1], cgp[:, 0]], 1) / cR[:, None] / np.hypot(cBp, cBt)[:, None]  # B_pol / |B|

    kind_cell = np.zeros(ntot, int)
    kind_cell[R[bnd]] = face_kind[bnd]

    # Jacobian sparsity: cells within two face hops, through real cells only
    nb = [set() for _ in range(nreal)]
    for a_, b_ in zip(L, R):
        if b_ < nreal:
            nb[a_].add(b_)
            nb[b_].add(a_)
    rows = []
    # Real cells in the vertex stencils of the faces of each cell
    stencil_cells = [set() for _ in range(nreal)]
    for f in range(nf):
        st = [int(c) for c in scell[sptr[f]:sptr[f + 1]] if c < nreal]
        stencil_cells[L[f]].update(st)
        if R[f] < nreal:
            stencil_cells[R[f]].update(st)
    at_vertex = [[] for _ in range(len(vx))]
    for i, c in enumerate(cells):
        for v in c:
            at_vertex[v].append(i)
    for i in range(nreal):
        s = {i} | nb[i]
        for j in nb[i]:
            s |= nb[j]
        for v in cells[i]:      # the skew correction couples cells that share a vertex
            s.update(at_vertex[v])
        s.update(stencil_cells[i])
        for j in (prev[i], nxt[i]):
            if j >= 0:
                s.add(j)
                for k in (prev[j], nxt[j]):
                    if k >= 0:
                        s.add(k)
        rows.append({int(j) for j in s})
    # The wider stencils at wall vertices are one-sided: make the pattern symmetric
    for i in range(nreal):
        for j in rows[i]:
            rows[j].add(i)
    rows = [sorted(r) for r in rows]
    for g_ in range(nreal, ntot):
        rows.append([g_])
    rowptr = np.concatenate([[0], np.cumsum([len(r) for r in rows])])
    cols = np.concatenate(rows)

    nvmax = max(len(c) for c in cells)
    cell_vx = -np.ones((nreal, nvmax), int)
    for i, c in enumerate(cells):
        cell_vx[i, :len(c)] = c

    return dict(
        nreal=nreal, nghost=nghost, nfaces=nf, nclosed_tubes=nclosed,
        R=xc[:, 0], Z=xc[:, 1],
        vol=cellfield(vol, 1.0), area=cellfield(area, 1.0),
        bR=cellfield(cbvec[:, 0], 0.0), bZ=cellfield(cbvec[:, 1], 0.0),
        B=cellfield(np.hypot(cBp, cBt), Bf[bnd]), Bp=cellfield(cBp, Bp[bnd]), Bt=cellfield(cBt, Bt[bnd]),
        psi=cellfield(cpsi, 0.5 * (vpsi[v1] + vpsi[v2])[bnd]),
        tube=cellfield(tube_id, -1).astype(int), prev=cellfield(prev, -1).astype(int),
        next=cellfield(nxt, -1).astype(int), kind=kind_cell,
        face_L=L, face_R=R, face_flux=flux, face_B=Bf, face_Bp=Bp, face_area=farea, face_dn=dn,
        face_length=length, face_Rmid=mid[:, 0], face_Zmid=mid[:, 1],
        face_kpar=kpar, face_tperp=tperp, face_tiso=tiso, face_kind=face_kind,
        face_drift=fdrift, face_curl=fcurl,
        face_vperp=vperp, face_viso=viso, face_sptr=sptr, face_scell=scell, face_sw=sw,
        vertex_kind=vertex_kind,
        face_label=np.zeros(nf, int) if face_label is None else np.asarray(face_label, int),
        face_side_L=side_L, face_side_R=side_R, face_up=up,
        face_v1=v1, face_v2=v2, face_nR=nrm[:, 0], face_nZ=nrm[:, 1],
        jac_rowptr=rowptr, jac_cols=cols,
        vertex_R=vx, vertex_Z=vy, cell_vertices=cell_vx,
    )


def _vertex_stencils(P, xc, cells, v1, v2, L, R, face_kind, nreal):
    """Difference of the two vertex values of every face as a weighted sum over cells.

    The value at a vertex is interpolated from cells near it:
      core vertex       the ghost cells of the core faces that meet there
                        (they hold the boundary value), weighted by distance
      other vertices    a weighted least-squares linear fit through the real
                        cells that share the vertex, which reproduces a linear
                        function exactly. If those cells do not determine the
                        fit (fewer than three, or nearly on a line, as happens
                        at the wall) their face neighbours are added.
                        Wall ghost cells are not used: not every field keeps
                        them up to date.
    Returns CSR arrays (ptr, cell, weight) with one row per face, and the
    kind of each vertex (0 interior, 1 core, 2 wall).
    """
    nv = len(P)
    around = [[] for _ in range(nv)]
    for i, c in enumerate(cells):
        for v in c:
            around[v].append(i)
    neighbours = [set() for _ in range(nreal)]
    for a_, b_ in zip(L, R):
        if b_ < nreal:
            neighbours[a_].add(b_)
            neighbours[b_].add(a_)
    core_ghosts = [[] for _ in range(nv)]
    vkind = np.zeros(nv, int)
    for f in np.flatnonzero(face_kind != KIND_INTERIOR):
        for v in (v1[f], v2[f]):
            if face_kind[f] == KIND_CORE:
                core_ghosts[v].append(R[f])
                vkind[v] = KIND_CORE
            elif vkind[v] != KIND_CORE:
                vkind[v] = KIND_WALL

    def fit(ids, x0):
        """weights of a weighted least-squares linear fit evaluated at x0, or None"""
        r = xc[ids] - x0
        h = np.hypot(*r.T).max()
        A = np.column_stack([np.ones(len(ids)), r / h])
        w = 1.0 / np.maximum((r**2).sum(1) / h**2, 1e-12)
        M = A.T @ (w[:, None] * A)
        if len(ids) < 3 or np.linalg.cond(M) > 1e6:
            return None
        wts = np.linalg.solve(M, A.T * w)[0]
        return wts if np.abs(wts).max() < 3.0 else None

    weights = []
    for v in range(nv):
        if vkind[v] == KIND_CORE:
            ids = np.array(core_ghosts[v])
            dist = np.maximum(np.hypot(*(xc[ids] - P[v]).T), 1e-300)
            weights.append((ids, (1 / dist) / (1 / dist).sum()))
            continue
        ids = sorted(around[v])
        if not ids:
            weights.append((np.zeros(0, int), np.zeros(0)))
            continue
        wts = fit(np.array(ids), P[v])
        if wts is None:
            wider = set(ids)
            for i in ids:
                wider |= neighbours[i]
            wide = np.array(sorted(wider))
            wts = fit(wide, P[v])
            if wts is not None:
                ids = wide
        if wts is None:      # last resort: inverse distance
            dist = np.maximum(np.hypot(*(xc[np.array(ids)] - P[v]).T), 1e-300)
            wts = (1 / dist) / (1 / dist).sum()
            vkind[v] += 10   # marks a vertex without a linear fit
        weights.append((np.array(ids), wts))
    ptr, cell, sw = [0], [], []
    for f in range(len(v1)):
        # Wall faces have a row too: no gradient is taken across them (their
        # transmissibilities are zero), but the E x B flow into the wall needs
        # the potential difference along the face.
        row = {}
        for v, sign in ((v2[f], 1.0), (v1[f], -1.0)):
            for i, w in zip(*weights[v]):
                row[int(i)] = row.get(int(i), 0.0) + sign * w
        for i in sorted(row):
            if abs(row[i]) > 1e-14:
                cell.append(i)
                sw.append(row[i])
        ptr.append(len(cell))
    return np.array(ptr), np.array(cell, int), np.array(sw), vkind


def check(g, verbose=True, strict=True):
    """Consistency checks. Returns a dict of numbers; raises if the grid is broken."""
    nreal, L, R = g["nreal"], g["face_L"], g["face_R"]
    out = {}
    # discrete div B = 0
    div = np.zeros(g["nreal"] + g["nghost"])
    np.add.at(div, L, g["face_flux"])
    np.add.at(div, R, -g["face_flux"])
    scale = np.zeros_like(div)
    np.add.at(scale, L, np.abs(g["face_flux"]))
    np.add.at(scale, R, np.abs(g["face_flux"]))
    out["max |sum of face fluxes| / sum |flux| per real cell"] = float(np.max(np.abs(div[:nreal]) / np.maximum(scale[:nreal], 1e-300)))
    # closed surface: sum of outward area vectors vanishes in the (R, Z) plane
    ax = np.zeros_like(div)
    az = np.zeros_like(div)
    length = g["face_area"] / (0.5 * (g["vertex_R"][g["face_v1"]] + g["vertex_R"][g["face_v2"]]))
    for arr, n in ((ax, g["face_nR"]), (az, g["face_nZ"])):
        np.add.at(arr, L, n * length)
        np.add.at(arr, R, -n * length)
    per = np.zeros_like(div)
    np.add.at(per, L, length)
    np.add.at(per, R, length)
    out["max |sum of face normals * length| / perimeter per real cell"] = float(np.max(np.hypot(ax, az)[:nreal] / per[:nreal]))
    out["min face normal distance [m]"] = float(g["face_dn"].min())
    out["faces with non-positive normal distance"] = int((g["face_dn"] <= 0).sum())
    out["min cell volume [m^3/rad]"] = float(g["vol"][:nreal].min())
    out["total volume [m^3]"] = float(2 * np.pi * g["vol"][:nreal].sum())
    out["cells without a tube neighbour"] = int(((g["prev"][:nreal] < 0) & (g["next"][:nreal] < 0)).sum())
    out["closed flux tubes"] = int(g.get("nclosed_tubes", -1))
    tube_faces = (g["face_side_L"] != 0)
    out["faces joining tube neighbours"] = int(tube_faces.sum())
    out["median |b_pol.n| on tube faces / on other interior faces"] = (
        float(np.median(np.abs(g["face_flux"] / (g["face_area"] * g["face_Bp"]))[tube_faces])),
        float(np.median(np.abs(g["face_flux"] / (g["face_area"] * g["face_Bp"]))[(~tube_faces) & (g["face_kind"] == 0)])))
    # Parallel transport must not leak through faces that carry cross-field
    # transport: parallel heat conduction is ~1e10 times stronger at 1 keV
    cross = (g["face_kind"] == KIND_INTERIOR) & (g["face_side_L"] == 0) & (g["face_tperp"] > 0)
    leak = g["face_kpar"][cross] / g["face_tperp"][cross] if cross.any() else np.zeros(1)
    out["parallel / cross-field transmissibility on flux-surface faces (max; must be < 1e-12)"] = float(leak.max())
    out["boundary faces core / wall"] = (int((g["face_kind"] == KIND_CORE).sum()), int((g["face_kind"] == KIND_WALL).sum()))
    out["wall faces with an upstream cell"] = int(((g["face_kind"] == KIND_WALL) & (g["face_up"] >= 0)).sum())
    out["wall area crossed by the field, sum |flux|/B [m^2]"] = float(2 * np.pi * (np.abs(g["face_flux"]) / g["face_B"])[g["face_kind"] == KIND_WALL].sum())
    if "jac_rowptr" in g:
        out["Jacobian non-zeros per real cell (mean, max)"] = (float(np.diff(g["jac_rowptr"])[:nreal].mean()), int(np.diff(g["jac_rowptr"]).max()))
    if verbose:
        for k, v in out.items():
            print("   %-66s %s" % (k, v))
    assert out["faces with non-positive normal distance"] == 0
    if strict:
        assert out["max |sum of face fluxes| / sum |flux| per real cell"] < 1e-10
        assert leak.max() < 1e-12, "parallel transport leaks through flux-surface faces"
    return out


def layout(g, blocks=48, order="angle"):
    """Storage order of the cells in the grid file.

    Real cells are sorted by poloidal angle (order="angle") or kept flux tube
    by flux tube (order="tube"). Each real cell is followed directly by its
    ghost cells, so that a cell and its boundary values are always on the
    same processor. Inert padding
    cells are inserted so that the total is a multiple of `blocks` and no
    such group straddles a multiple of total/blocks: BOUT++ splits the list
    into equal contiguous pieces, so any processor count that divides
    `blocks` then works.

    Returns (order, storage_of_logical): order[s] is the logical cell stored
    at position s, or -1 for padding.
    """
    nreal, ntot = g["nreal"], g["nreal"] + g["nghost"]
    groups = [[c] for c in range(nreal)]
    for L, R in zip(g["face_L"], g["face_R"]):
        if R >= nreal:
            groups[L].append(R)
    if order == "angle":
        # Sort by poloidal angle about the middle of the core boundary, so that
        # each contiguous piece is a poloidal sector with short interfaces.
        # (Flux-tube order gives long thin pieces whose halo is as large as
        # the piece itself.)
        core = g["face_kind"] == KIND_CORE
        if core.any():
            R0, Z0 = g["face_Rmid"][core].mean(), g["face_Zmid"][core].mean()
        else:
            R0, Z0 = g["R"][:nreal].mean(), g["Z"][:nreal].mean()
        angle = np.arctan2(g["Z"][:nreal] - Z0, g["R"][:nreal] - R0)
        groups = [groups[c] for c in np.argsort(angle, kind="stable")]
    elif order != "tube":
        raise ValueError("unknown cell order %r" % order)
    longest = max(len(grp) for grp in groups)
    N = -(-ntot // blocks) * blocks
    while True:
        b = N // blocks
        order = []
        if b >= longest:
            for grp in groups:
                used = len(order) % b
                if used + len(grp) > b:
                    order.extend([-1] * (b - used))
                order.extend(grp)
            if len(order) <= N:
                order.extend([-1] * (N - len(order)))
                break
        N += blocks
    order = np.array(order)
    storage_of_logical = np.empty(ntot, int)
    storage_of_logical[order[order >= 0]] = np.flatnonzero(order >= 0)
    return order, storage_of_logical


# Value given to padding cells
_PAD = {"vol": 1.0, "area": 1.0, "B": 1.0, "kind": KIND_PADDING, "prev": -1, "next": -1, "tube": -1}
# Arrays whose values are cell indices
_CELL_INDEX = ("prev", "next", "face_L", "face_R", "face_up", "face_scell")


def write(g, path, extra=None, blocks=48, order="angle"):
    """Write a BOUT++ grid file. `extra` holds further per-cell arrays (logical order)."""
    ntot = g["nreal"] + g["nghost"]
    order, s_of_l = layout(g, blocks, order)
    N = len(order)
    real = order >= 0

    def cells(name, a):
        """Per-cell array from logical to storage order"""
        a = np.asarray(a)
        out = np.full(N, _PAD.get(name, 0), dtype=a.dtype)
        out[real] = a[order[real]]
        if name in _CELL_INDEX:
            out[real] = np.where(out[real] >= 0, s_of_l[np.maximum(out[real], 0)], -1)
        return out

    def remap(a):
        a = np.asarray(a)
        return np.where(a >= 0, s_of_l[np.maximum(a, 0)], -1)

    # Jacobian pattern rows in storage order
    rows = []
    for s in range(N):
        c = order[s]
        if c < 0:
            rows.append([s])
        else:
            rows.append(sorted(s_of_l[g["jac_cols"][g["jac_rowptr"][c]:g["jac_rowptr"][c + 1]]].tolist()))
    rowptr = np.concatenate([[0], np.cumsum([len(r) for r in rows])])
    cols = np.concatenate(rows)

    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", 1)
        d.createDimension("y", N)
        d.createDimension("us_cell", N)
        d.createDimension("us_logical", ntot)
        d.createDimension("us_face", g["nfaces"])
        d.createDimension("us_jac", len(cols))
        d.createDimension("us_jacrow", N + 1)
        d.createDimension("us_vertex", len(g["vertex_R"]))
        d.createDimension("us_real", g["nreal"])
        d.createDimension("us_nv", g["cell_vertices"].shape[1])

        def scalar(name, val):
            v = d.createVariable(name, "i4")
            v[...] = val

        scalar("nx", 1)
        scalar("ny", N)
        scalar("ixseps1", -1)
        scalar("ixseps2", -1)
        scalar("us_ncells", g["nreal"])
        scalar("us_nghost", g["nghost"])
        scalar("us_npad", N - ntot)
        scalar("us_nfaces", g["nfaces"])
        scalar("us_blocks", blocks)

        # Placeholder metric: Hermes-3 replaces J and Bxy when it loads the cells
        for name, val in (("dx", 1.0), ("dy", 1.0), ("g11", 1.0), ("g22", 1.0), ("g33", 1.0), ("g12", 0.0),
                          ("g13", 0.0), ("g23", 0.0), ("J", 1.0), ("Bxy", 1.0)):
            v = d.createVariable(name, "f8", ("x", "y"))
            v[0, :] = val
        v = d.createVariable("us_logical", "i4", ("us_logical",))
        v[:] = s_of_l
        for name in ("R", "Z", "vol", "B", "Bp", "Bt", "psi", "area", "bR", "bZ"):
            v = d.createVariable("us_" + name, "f8", ("us_cell",))
            v[:] = cells(name, g[name])
        for name in ("tube", "prev", "next", "kind"):
            v = d.createVariable("us_" + name, "i4", ("us_cell",))
            v[:] = cells(name, g[name])
        for name, val in (extra or {}).items():
            v = d.createVariable("us_" + name, "f8", ("us_cell",))
            v[:] = cells(name, val)
        for name in ("flux", "B", "Bp", "area", "length", "dn", "kpar", "tperp", "tiso", "nR", "nZ", "Rmid", "Zmid",
                     "drift", "curl"):
            if "face_" + name not in g:
                continue
            v = d.createVariable("us_face_" + name, "f8", ("us_face",))
            v[:] = g["face_" + name]
        for name in ("L", "R", "kind", "label", "side_L", "side_R", "up", "v1", "v2"):
            v = d.createVariable("us_face_" + name, "i4", ("us_face",))
            v[:] = remap(g["face_" + name]) if "face_" + name in _CELL_INDEX else g["face_" + name]
        if "face_sptr" in g:
            d.createDimension("us_faceptr", g["nfaces"] + 1)
            d.createDimension("us_stencil", len(g["face_scell"]))
            for name in ("vperp", "viso"):
                v = d.createVariable("us_face_" + name, "f8", ("us_face",))
                v[:] = g["face_" + name]
            v = d.createVariable("us_face_sptr", "i4", ("us_faceptr",))
            v[:] = g["face_sptr"]
            v = d.createVariable("us_face_scell", "i4", ("us_stencil",))
            v[:] = remap(g["face_scell"])
            v = d.createVariable("us_face_sw", "f8", ("us_stencil",))
            v[:] = g["face_sw"]
            v = d.createVariable("us_vertex_kind", "i4", ("us_vertex",))
            v[:] = g["vertex_kind"]
        v = d.createVariable("us_jac_rowptr", "i4", ("us_jacrow",))
        v[:] = rowptr
        v = d.createVariable("us_jac_cols", "i4", ("us_jac",))
        v[:] = cols
        for name in ("R", "Z"):
            v = d.createVariable("us_vertex_" + name, "f8", ("us_vertex",))
            v[:] = g["vertex_" + name]
        v = d.createVariable("us_cell_vertices", "i4", ("us_real", "us_nv"))
        v[:] = g["cell_vertices"]
    return N


def read(path):
    """Read the unstructured part of a grid file, back in logical order:
    real cells first (flux tube by flux tube), then the ghost cells.
    `storage_index` gives the position of each logical cell in the file and
    in Hermes-3 output."""
    with netCDF4.Dataset(path) as d:
        s_of_l = np.asarray(d["us_logical"][:])
        N = len(d.dimensions["us_cell"])
        l_of_s = -np.ones(N, int)
        l_of_s[s_of_l] = np.arange(len(s_of_l))
        g = {}
        for k, v in d.variables.items():
            if not k.startswith("us_") or k == "us_logical" or k.startswith("us_jac"):
                continue
            a = np.asarray(v[...])
            name = k[3:]
            if v.dimensions == ("us_cell",):
                a = a[s_of_l]
            if name in _CELL_INDEX:
                a = np.where(a >= 0, l_of_s[np.maximum(a, 0)], -1)
            g[name] = a
    g["nreal"], g["nghost"], g["nfaces"] = int(g.pop("ncells")), int(g["nghost"]), int(g["nfaces"])
    g["storage_index"] = s_of_l
    g["nstorage"] = N
    return g
