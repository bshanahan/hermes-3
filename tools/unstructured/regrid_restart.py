"""Copy the state of a flat-cell run onto another flat-cell grid.

usage: python regrid_restart.py <source dir> <source grid> <target dir> <target grid>

Every evolving variable in the restart files of <target dir> is overwritten
with values of the source run: a real cell takes the value of the source cell
that, among the 16 nearest, lies on the closest flux surface; a ghost cell
that of the nearest ghost cell of the same kind (by the mid-point of its
boundary face). Nothing is interpolated across flux surfaces. The target restart
files must exist already, with the processor count that will be used: make
them with a very short run from the same input (nout = 1 and a tiny timestep).

This gives a starting state for a run on a refined or coarsened grid of the
same geometry. It is not conservative.
"""
import glob
import sys

import netCDF4
import numpy as np
from scipy.spatial import cKDTree

import flatgrid


def restart_files(path):
    files = sorted(glob.glob(path + "/BOUT.restart.*.nc"), key=lambda f: int(f.split(".")[-2]))
    if not files:
        raise SystemExit("no restart files in " + path)
    return files


def read_state(path, g):
    """{name: values in logical order} for every evolving variable of a run"""
    files = [netCDF4.Dataset(f) for f in restart_files(path)]
    myg = int(files[0]["MYG"][...])
    out = {}
    for name, v in files[0].variables.items():
        if v.dimensions != ("x", "y", "z"):
            continue
        a = np.concatenate([np.asarray(f[name][0, myg:f[name].shape[1] - myg, 0]) for f in files])
        if len(a) != g["nstorage"]:
            raise SystemExit("%s: %d values in the restart files, grid has %d cells" % (name, len(a), g["nstorage"]))
        out[name] = a[g["storage_index"]]
    for f in files:
        f.close()
    return out


def main(src_dir, src_grid, dst_dir, dst_grid):
    gs, gd = flatgrid.read(src_grid), flatgrid.read(dst_grid)
    state = read_state(src_dir, gs)
    ns, nd = gs["nreal"], gd["nreal"]

    # Nearest source cell of every target cell, kind by kind
    source_of = np.empty(nd + gd["nghost"], int)
    for kind in (flatgrid.KIND_INTERIOR, flatgrid.KIND_CORE, flatgrid.KIND_WALL):
        s = np.flatnonzero(gs["kind"] == kind)
        d = np.flatnonzero(gd["kind"] == kind)
        if len(d) == 0:
            continue
        if len(s) == 0:
            raise SystemExit("the source grid has no cells of kind %d" % kind)
        dist, j = cKDTree(np.stack([gs["R"][s], gs["Z"][s]], 1)).query(np.stack([gd["R"][d], gd["Z"][d]], 1))
        source_of[d] = s[j]
        print("kind %d: %d target cells, distance to the source cell median %.1f mm, max %.1f mm" % (kind, len(d), 1e3 * np.median(dist), 1e3 * dist.max()))

    # Real cells: of the nearest source cells, the one on the closest flux
    # surface. Interpolating across flux surfaces would give the potential and
    # the pressures a gradient along the field, and with it very large
    # parallel currents and heat flows.
    real_s = np.flatnonzero(gs["kind"] == flatgrid.KIND_INTERIOR)
    real_d = np.flatnonzero(gd["kind"] == flatgrid.KIND_INTERIOR)
    k = min(16, len(real_s))
    dist, j = cKDTree(np.stack([gs["R"][real_s], gs["Z"][real_s]], 1)).query(
        np.stack([gd["R"][real_d], gd["Z"][real_d]], 1), k=k)
    dpsi = np.abs(gs["psi"][real_s][j] - gd["psi"][real_d][:, None])
    # distance only breaks ties between cells of the same surface
    score = dpsi / max(np.ptp(gs["psi"][real_s]), 1e-300) + 1e-6 * dist / dist.max()
    best = np.argmin(score, axis=1)
    rows = np.arange(len(real_d))
    source_of[real_d] = real_s[j[rows, best]]
    print("real cells: psi mismatch to the source cell median %.1e, max %.1e of the psi range; distance median %.1f mm, max %.1f mm" % (
        np.median(dpsi[rows, best]) / np.ptp(gs["psi"][real_s]), dpsi[rows, best].max() / np.ptp(gs["psi"][real_s]),
        1e3 * np.median(dist[rows, best]), 1e3 * dist[rows, best].max()))

    def target_values(values):
        return values[source_of]

    files = restart_files(dst_dir)
    nlocal = gd["nstorage"] // len(files)
    logical_of_storage = -np.ones(gd["nstorage"], int)
    logical_of_storage[gd["storage_index"]] = np.arange(len(gd["storage_index"]))
    done = set()
    target = {name: target_values(values) for name, values in state.items()}
    for rank, path in enumerate(files):
        with netCDF4.Dataset(path, "a") as f:
            myg = int(f["MYG"][...])
            logical = logical_of_storage[rank * nlocal:(rank + 1) * nlocal]
            used = logical >= 0
            for name, source_values in state.items():
                values = target[name]
                if name not in f.variables:
                    v = f.createVariable(name, "f8", ("x", "y", "z"))
                    template = next(u for u in f.variables.values() if u.dimensions == ("x", "y", "z") and u is not v)
                    for a in template.ncattrs():
                        v.setncattr(a, template.getncattr(a))
                    v[:] = 0.0
                v = f[name]
                if v.shape[1] != nlocal + 2 * myg:
                    raise SystemExit("%s in %s has %d y points, expected %d" % (name, path, v.shape[1], nlocal + 2 * myg))
                a = np.asarray(v[0, :, 0])
                block = a[myg:myg + nlocal]
                block[used] = values[logical[used]]
                a[myg:myg + nlocal] = block
                v[0, :, 0] = a
                done.add(name)
    print("wrote %s into %d restart files of %s" % (", ".join(sorted(done)), len(files), dst_dir))


if __name__ == "__main__":
    main(*sys.argv[1:5])
