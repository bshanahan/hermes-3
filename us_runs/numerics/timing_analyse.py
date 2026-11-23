"""Approach to steady state against wall-clock time for the TCV-X21 timing runs.

usage (from us_runs/): python numerics/timing_analyse.py out.json <run dir>[@@grid] ...
       (a run dir without @@ is flat if its BOUT.inp names a flat grid, else structured)

For every output: plasma time, wall time, the distance to the final state
    d = max over n_i, p_i, p_e, n_n of ||F - F_final||_V / ||F_final||_V
(volume-weighted 2-norms over the real cells), and the relative rate of change
between outputs. Also: number of SNES steps and failures from BOUT.log.0.
"""
import json
import re
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "unstructured"))

FIELDS = ("Nd+", "Pd+", "Pe", "Nd")


def flat_series(run, grid):
    import usdata
    r = usdata.Run(run, grid)
    w = r.g["vol"][:r.n]
    data = {f: np.asarray(r(f))[:, :r.n] for f in FIELDS}
    for f in r.files:      # a second open of the same files in this process crashes HDF5
        f.close()
    with netCDF4.Dataset(run + "/BOUT.dmp.0.nc") as d:
        wall = np.asarray(d["wall_time"][:])
        t = np.asarray(d["t_array"][:]) / float(d["Omega_ci"][...])
    return t, wall, data, w, r.n


def struct_series(run):
    files = sorted(Path(run).glob("BOUT.dmp.*.nc"), key=lambda p: int(p.name.split(".")[2]))
    data = {f: [] for f in FIELDS}
    wts = []
    for p in files:
        with netCDF4.Dataset(str(p)) as d:
            nx, ny = d["MXSUB"][...], d["MYSUB"][...]
            mxg, myg = int(d["MXG"][...]), int(d["MYG"][...])
            pe_x = int(d["PE_XIND"][...])
            nxpe = int(d["NXPE"][...])
            xs = mxg if pe_x > 0 else 0
            xe = (mxg + nx) if pe_x < nxpe - 1 else (2 * mxg + nx)
            # real (non-boundary) cells only: drop x guard cells everywhere
            xs, xe = mxg, mxg + nx
            ys, ye = myg, myg + ny
            J = np.asarray(d["J"][:])[xs:xe, ys:ye]
            dx = np.asarray(d["dx"][:])[xs:xe, ys:ye]
            dy = np.asarray(d["dy"][:])[xs:xe, ys:ye]
            wts.append((J * dx * dy).ravel())
            for f in FIELDS:
                a = np.asarray(d[f][:])
                data[f].append(a[:, xs:xe, ys:ye, 0].reshape(a.shape[0], -1))
            if p.name == "BOUT.dmp.0.nc":
                wall = np.asarray(d["wall_time"][:])
                t = np.asarray(d["t_array"][:]) / float(d["Omega_ci"][...])
    data = {f: np.concatenate(v, 1) for f, v in data.items()}
    w = np.concatenate(wts)
    return t, wall, data, w, len(w)


def analyse(arg):
    run, grid = (arg.split("@@") + [None])[:2]
    inp = Path(run, "BOUT.inp").read_text()
    if grid is None:
        m = re.search(r'^file\s*=\s*"([^"]+)"', inp, re.M)
        grid = m.group(1) if m else None
        log = Path(run, "BOUT.log.0").read_text(errors="replace")
        m2 = re.search(r"Command line options for this run :.*mesh:file=(\S+)", log)
        if m2:
            grid = m2.group(1).strip('"')
    flat = "us_ncells" in netCDF4.Dataset(grid).variables
    t, wall, data, w, ncells = flat_series(run, grid) if flat else struct_series(run)
    dist = np.zeros(len(t))
    rate = np.full(len(t), np.nan)
    for f, a in data.items():
        nrm = np.sqrt(np.sum(w * a[-1]**2))
        dist = np.maximum(dist, np.sqrt(np.sum(w * (a - a[-1])**2, 1)) / nrm)
        ch = np.sqrt(np.sum(w * np.diff(a, axis=0)**2, 1)) / nrm / np.diff(t)
        rate[1:] = np.fmax(rate[1:], ch)
    log = Path(run, "BOUT.log.0").read_text(errors="replace")
    steps = len(re.findall(r"reason: \d", log))
    fails = len(re.findall(r"reason: -\d", log))
    ranks = int(re.search(r"Processor number: 0 of (\d+)", log).group(1))
    nl = [int(v) for v in re.findall(r"nl iter: (\d+), lin iter: \d+, reason: \d", log)]
    lin = [int(v) for v in re.findall(r"nl iter: \d+, lin iter: (\d+), reason: \d", log)]
    return {"run": run, "grid": grid, "flat": flat, "cells": int(ncells), "ranks": ranks,
            "unknowns": int(ncells) * 7, "steps": steps, "failures": fails,
            "newton_per_step": float(np.mean(nl)) if nl else None,
            "krylov_per_newton": float(np.sum(lin) / max(np.sum(nl), 1)) if nl else None,
            "t_ms": (1e3 * t).tolist(), "wall_s": wall.tolist(), "distance": dist.tolist(),
            "rate_per_s": rate.tolist()}


def time_to(res, level):
    """wall time [s] and plasma time [ms] after which the distance stays below `level`"""
    d = np.asarray(res["distance"])
    above = np.flatnonzero(d >= level)
    k = above[-1] + 1 if len(above) else 0
    if k >= len(d):
        return None, None
    return res["wall_s"][k], res["t_ms"][k]


if __name__ == "__main__":
    out = [analyse(a) for a in sys.argv[2:]]
    for r in out:
        w1, t1 = time_to(r, 1e-2)
        w3, t3 = time_to(r, 1e-3)
        print("%-34s %s cells %5d ranks %2d steps %5d fail %3d  Newton/step %.2f Krylov/Newton %.2f  total %.0f s | 1%%: %.0f s (%.0f ms) | 0.1%%: %.0f s (%.0f ms)"
              % (r["run"], "flat  " if r["flat"] else "struct", r["cells"], r["ranks"], r["steps"], r["failures"],
                 r["newton_per_step"] or 0, r["krylov_per_newton"] or 0, r["wall_s"][-1], w1 or np.nan, t1 or np.nan,
                 w3 or np.nan, t3 or np.nan))
    json.dump(out, open(sys.argv[1], "w"))
