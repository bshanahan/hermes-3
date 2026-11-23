"""Headline numbers of several TCV-X21 runs side by side, each on its own GOAT grid.

usage (from us_runs/): python tcvx21/scan.py <label>=<run dir>::<grid file> [...]
       (":: grid file" may be left out for the standard grid)

For a resolution scan: midplane separatrix values, flows, target peaks and the
distance of every profile to the forward-field measurements.
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np
import compare as c
import reference

twopi = 2*np.pi
c.SELECT = "tube"       # one cell per flux tube along each measurement line, at any resolution
STANDARD_GRID = c.GRID
rows = []
for arg in sys.argv[1:]:
    label, rest = arg.split("=", 1)
    run, grid = (rest.replace("@@", "::").split("::") + [STANDARD_GRID])[:2]
    c.GRID = grid
    d = c.load(run); g, n = d["g"], d["n"]
    out = {"cells": n, "t [ms]": 1e3*float(d["t"])}
    vals = reference.numbers(run)
    for k in ("midplane separatrix n_e [m^-3]", "midplane separatrix T_e [eV]", "midplane separatrix T_i [eV]", "ion flow to walls through the sheath [1/s]",
              "sheath power [W]", "ion content", "atom content"):
        out[k] = vals[k]
    # Midplane separatrix values from the flux label of each cell in the grid file
    # (the mean over its vertices). The label of the cell centre, which the
    # profiles use, moves off the flux surface as the cells get longer.
    omp = c.midplane_cells(g, n); omp = omp[np.argsort(d["grid_psin"][omp])]
    xg = c.rsep_of_psin(d["grid_psin"][omp])
    for k, a in (("midplane separatrix n_e [m^-3]", d["Ne"]), ("midplane separatrix T_e [eV]", d["Te"]), ("midplane separatrix T_i [eV]", d["Ti"])):
        out[k] = float(np.interp(0.0, xg, a[omp]))
    first = omp[np.searchsorted(xg, 0.0)]
    out["first open tube at the midplane: n_e [m^-3]"] = float(d["Ne"][first])
    out["first open tube at the midplane: T_e [eV]"] = float(d["Te"][first])
    fR, fZ = g["face_Rmid"], g["face_Zmid"]; fb = np.abs(d["fb"]); R = d["R"]
    for name, m in (("outer target", d["wall"] & (fZ < -0.749) & (fR > 0.69)), ("inner target", d["wall"] & (fR < 0.6245) & (fZ < -0.2))):
        out[name + ": ion flow [1/s]"] = float(twopi*(d["gam"]*fb)[m].sum())
        out[name + ": peak n_e [m^-3]"] = float(d["Ne"][R[m]].max())
        out[name + ": peak T_e [eV]"] = float(d["Te"][R[m]].max())
        out[name + ": peak J_sat [kA/m^2]"] = float(c.qe*d["gam"][m].max()/1e3)
    if d["r"].has("phi"):
        phi = d["r"]("phi", True)[-1]
        out["midplane separatrix phi [V]"] = float(np.interp(0.0, xg, phi[omp]))
        out["first open tube at the midplane: phi [V]"] = float(phi[first])
        out["lowest phi [V]"] = float(phi[:n].min())
    dj = [v for k, v in vals.items() if k.startswith("d ") and k.endswith("forward")]
    out["mean distance to the forward-field data"] = float(np.nanmean(dj))
    for k, v in vals.items():
        if k.startswith("d ") and k.endswith("forward"):
            out[k] = v
    rows.append((label, out))

keys = []
for _, o in rows:
    keys += [k for k in o if k not in keys]
w = max(18, max(len(lab) for lab, _ in rows) + 2)
print("%-46s" % "" + "".join(("%" + str(w) + "s") % lab for lab, _ in rows))
for k in keys:
    print("%-46s" % k + "".join(("%" + str(w) + ".4g") % o.get(k, np.nan) for _, o in rows))
