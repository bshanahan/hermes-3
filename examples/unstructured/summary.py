"""Key results of a flat-cell (unstructured) Hermes-3 tokamak run, as JSON.

usage: python summary.py <run directory> <grid file> [--json out.json] [--figure out.png]
                         [--compare reference.json] [--index -1] [--gamma-e 4.5] [--gamma-i 2.5]

Everything is in SI units, temperatures in eV. Printed and written:

  time_s, outputs      simulated time of the output used, number of outputs
  change_per_ms        relative change per millisecond between the last two outputs of the
                       electron content, the neutral content and the stored thermal energy
                       (all three small = steady)
  history              the same three totals at every output (written to the JSON file only)
  power_W              heating, sheath loss (electrons, ions), impurity radiation, neutral
                       energy to the walls, energy returned with recycled neutrals, and the
                       remainder (hydrogen radiation, charge exchange and the rate of change)
  particles_per_s      ion flow to the walls
  core_boundary        volume average over the cells that touch the core boundary
  outer_midplane       separatrix values and the profile along the outboard midplane
                       (psi in Wb/rad, zero on the separatrix, negative inside)
  targets              per divertor target: power, ion flow, and at the wall face with the
                       largest heat flux its position, the heat flux along the field, the heat
                       flux on the surface and the plasma state in front of it

The heat flux on a wall face is rebuilt from the sheath-edge values that the sheath
boundary stores in the ghost cell of the face,
    q_par = gamma_e T_e n_e v_e + sum over ions of gamma_i T_i n_i v_i      (along B)
    q_surface = (q_par + 13.6 eV * ion flux) * |b.n|                        (on the surface)
with gamma_e and gamma_i as in [sheath_boundary_unstructured] of the input. The sum over
all faces is checked against the Ee_sheath and Ei_sheath diagnostics
(sheath_check_faces_over_diagnostic should be 1). Radiation and neutral loads on the
surface are not included.

Targets are found from the grid alone: an X-point is where the poloidal field of the
cells is smallest in the lower (upper) half; the wall faces below (above) it are the
lower (upper) divertor, split into inner and outer at its major radius. The rest of the
wall is the main chamber.

--compare prints the main quantities next to those of a reference summary (for
example sparc-lmode-state/summary.json) with their relative differences.

Needs numpy and netCDF4 (and matplotlib for --figure).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools", "unstructured"))
import usdata  # noqa: E402

QE = 1.602176634e-19
TWOPI = 2 * np.pi


def find_species(run):
    """Ions and neutrals with an evolved density, from the output variable names."""
    names = [v[1:] for v in run.d.variables
             if v.startswith("N") and not v.startswith("NV") and run.d[v].ndim == 4 and v != "Ne"]
    ions = [s for s in names if s.endswith("+")]
    neutrals = [s for s in names if not s.endswith("+")]
    return ions, neutrals


def xpoints(g, n):
    """(R, Z) of the lower and upper X-point, or None: where B_pol of the cells is smallest."""
    Z, R, Bp = g["Z"][:n], g["R"][:n], g["Bp"][:n]
    zmid = 0.5 * (Z.min() + Z.max())
    out = []
    for half in (Z < zmid, Z >= zmid):
        idx = np.flatnonzero(half)
        c = idx[np.argmin(Bp[idx])]
        out.append((float(R[c]), float(Z[c])) if Bp[c] < 0.1 * np.median(Bp) else None)
    return out


def summarise(path, gridfile, index=-1, gamma_e=4.5, gamma_i=2.5):
    r = usdata.Run(path, gridfile)
    g, n, it = r.g, r.n, index
    ions, neutrals = find_species(r)
    Ni = {s: r("N" + s, True) for s in ions}
    Vi = {s: r("V" + s, True)[it] for s in ions}
    Ti = {s: r("T" + s, True)[it] for s in ions}
    Ne_all = sum(Ni.values())
    Ne = Ne_all[it]
    Te = r("Te", True)[it]
    Ve = r("Ve", True)[it]
    Pe_all = r("Pe", True)
    Pi_all = sum(r("P" + s, True) for s in ions)
    Nn_all = sum(r("N" + s, True) for s in neutrals) if neutrals else np.zeros_like(Ne_all)
    Tiavg = sum(Ni[s][it] * Ti[s] for s in ions) / np.maximum(Ne, 1.0)

    out = {"run": os.path.abspath(path), "grid": os.path.abspath(gridfile), "cells": int(n),
           "ions": ions, "neutrals": neutrals,
           "time_s": float(r.t[it]), "outputs": int(len(r.t))}

    # How steady: relative change per ms between the last two outputs
    content = r.volume_integral(Ne_all)
    neutral_content = r.volume_integral(Nn_all)
    energy = r.volume_integral(1.5 * (Pe_all + Pi_all))
    if len(r.t) > 1:
        dt_ms = 1e3 * (r.t[it] - r.t[it - 1])
        out["change_per_ms"] = {
            "interval_ms": float(dt_ms),
            "electron_content": float((content[it] / content[it - 1] - 1) / dt_ms),
            "neutral_content": float((neutral_content[it] / neutral_content[it - 1] - 1) / dt_ms)
            if neutral_content[it - 1] > 0 else 0.0,
            "stored_energy": float((energy[it] / energy[it - 1] - 1) / dt_ms)}
    out["electron_content"] = float(content[it])
    out["neutral_content"] = float(neutral_content[it])
    out["stored_energy_J"] = float(energy[it])
    out["history"] = {"time_s": [float(v) for v in r.t], "electron_content": content.tolist(),
                      "neutral_content": neutral_content.tolist(), "stored_energy_J": energy.tolist()}

    # Power balance
    def integral(name):
        return float(r.volume_integral(r(name, True)[it])) if r.has(name) else 0.0

    heating = sum(1.5 * integral("P%s_src" % s) for s in ions + ["e"])
    sheath_e, sheath_i = -integral("Ee_sheath"), -integral("Ei_sheath")
    rad_names = [v for v in r.d.variables if v.startswith("Rfixed_fraction") and r.d[v].ndim == 4]
    radiation = {v[1:]: abs(integral(v)) for v in rad_names}
    neutral_wall = -sum(integral("E%s_target_refl" % s) for s in neutrals)
    recycled = sum(integral("E%s_target_recycle" % s) for s in neutrals)
    out["power_W"] = {
        "heating": heating, "sheath_electrons": sheath_e, "sheath_ions": sheath_i,
        "impurity_radiation": radiation, "neutrals_to_walls": neutral_wall,
        "returned_with_recycled_neutrals": recycled,
        "remainder": heating - sheath_e - sheath_i - sum(radiation.values()) - neutral_wall + recycled}

    # Wall faces: sheath-edge values are in the ghost cell (index R) of each face
    L, R, kind = g["face_L"], g["face_R"], g["face_kind"]
    wall = kind == 2
    fb = g["face_flux"] / g["face_B"]              # area perpendicular to B per radian, signed
    sign = np.where(fb >= 0, 1.0, -1.0)
    flux_i = {s: np.where(wall, Ni[s][it][R] * np.maximum(sign * Vi[s][R], 0), 0.0) for s in ions}
    gam = sum(flux_i.values())                     # ions per m^2 per s along B
    q_par = np.where(wall, gamma_e * QE * Te[R] * Ne[R] * np.maximum(sign * Ve[R], 0), 0.0) \
        + sum(gamma_i * QE * Ti[s][R] * flux_i[s] for s in ions)
    bn = np.abs(fb) / g["face_area"]               # |b.n|
    q_surface = (q_par + 13.6 * QE * gam) * bn
    face_power = TWOPI * q_par * np.abs(fb)        # W per face
    face_flow = TWOPI * gam * np.abs(fb)           # ions/s per face
    out["particles_per_s"] = {"ions_to_walls": float(face_flow.sum())}
    out["sheath_check_faces_over_diagnostic"] = float(face_power.sum() / max(sheath_e + sheath_i, 1e-30))

    # Core boundary
    core_cells = np.unique(L[kind == 1])
    if len(core_cells):
        w = g["vol"][core_cells]

        def mean(a):
            return float((a[core_cells] * w).sum() / w.sum())

        out["core_boundary"] = {"psi_Wb_per_rad": mean(g["psi"]), "n_e": mean(Ne), "T_e": mean(Te),
                                "T_i": mean(Tiavg), "pressure_Pa": mean(Pe_all[it] + Pi_all[it])}

    # Outboard midplane: the cell of largest major radius in every flux tube whose
    # outermost point lies on the outboard side near the height of the widest tube
    tube = g["tube"][:n]
    Rc, Zc, psi = g["R"][:n], g["Z"][:n], g["psi"][:n]
    order = np.lexsort((Rc, tube))
    last = np.flatnonzero(np.r_[tube[order][1:] != tube[order][:-1], True])
    outer = order[last]                            # per tube: the cell with the largest R
    z0 = Zc[outer[np.argmax(Rc[outer])]]
    height = Zc.max() - Zc.min()
    outer = outer[(np.abs(Zc[outer] - z0) < 0.05 * height) & (Rc[outer] > 0.5 * (Rc.min() + Rc.max()))]
    outer = outer[np.argsort(psi[outer])]
    out["outer_midplane"] = {
        "separatrix": {k: float(np.interp(0.0, psi[outer], a[outer]))
                       for k, a in (("R", Rc), ("n_e", Ne), ("T_e", Te), ("T_i", Tiavg))},
        "profile": {"psi_Wb_per_rad": psi[outer].tolist(), "R": Rc[outer].tolist(), "Z": Zc[outer].tolist(),
                    "n_e": Ne[outer].tolist(), "T_e": Te[outer].tolist(), "T_i": Tiavg[outer].tolist()}}

    # Targets
    fR, fZ = g["face_Rmid"], g["face_Zmid"]
    lower, upper = xpoints(g, n)
    regions = {}
    assigned = np.zeros(len(L), bool)
    for name, xp, below in (("lower", lower, True), ("upper", upper, False)):
        if xp is None:
            continue
        side = wall & ((fZ < xp[1]) if below else (fZ > xp[1]))
        regions[name + "_inner"] = side & (fR <= xp[0])
        regions[name + "_outer"] = side & (fR > xp[0])
        assigned |= side
    regions["main_chamber"] = wall & ~assigned
    Nn = Nn_all[it]
    out["x_points"] = {"lower": lower, "upper": upper}
    out["targets"] = {}
    for name, sel in regions.items():
        f = np.flatnonzero(sel)
        if len(f) == 0:
            continue
        top = f[np.argmax(q_surface[f])]
        out["targets"][name] = {
            "power_W": float(face_power[f].sum()), "ion_flow_per_s": float(face_flow[f].sum()),
            "wall_faces": int(len(f)), "max_T_e": float(Te[L[f]].max()),
            "peak": {"R": float(fR[top]), "Z": float(fZ[top]),
                     "q_surface_W_per_m2": float(q_surface[top]), "q_parallel_W_per_m2": float(q_par[top]),
                     "T_e": float(Te[R[top]]), "n_e": float(Ne[R[top]]),
                     "neutral_density": float(Nn[L[top]])}}
    fields = dict(Ne=Ne[:n], Te=Te[:n], Ti=Tiavg[:n], Nn=Nn[:n],
                  rad=sum(np.abs(r(v)[it]) for v in rad_names) if rad_names else None)
    for f in r.files:   # HDF5 crashes if a process opens the same files again without closing them
        f.close()
    return out, g, fields


def figure(path, g, fields, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    panels = [(fields["Ne"], "electron density [m$^{-3}$]", "viridis", 1e-3),
              (fields["Te"], "electron temperature [eV]", "plasma", None),
              (fields["Nn"], "neutral density [m$^{-3}$]", "viridis", 1e-5)]
    if fields["rad"] is not None and fields["rad"].max() > 0:
        panels.append((fields["rad"], "impurity radiation [W m$^{-3}$]", "inferno", 1e-4))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 8.5))
    for ax, (val, label, cmap, floor) in zip(axes, panels):
        vmax = float(val.max())
        vmin = vmax * floor if floor else 0.5
        pc = usdata.plot_cells(ax, g, np.maximum(val, vmin), edgecolors="none", cmap=cmap,
                               norm=LogNorm(vmin, max(vmax, 2 * vmin)))
        fig.colorbar(pc, ax=ax, shrink=0.8)
        ax.set_title(label)
        ax.set_xlabel("R [m]")
    axes[0].set_ylabel("Z [m]")
    fig.suptitle(title)
    fig.savefig(path, dpi=110, bbox_inches="tight")


# Quantities compared by --compare: (label, path into the summary)
COMPARE = [
    ("separatrix n_e (outer midplane)", ("outer_midplane", "separatrix", "n_e")),
    ("separatrix T_e (outer midplane)", ("outer_midplane", "separatrix", "T_e")),
    ("separatrix T_i (outer midplane)", ("outer_midplane", "separatrix", "T_i")),
    ("core boundary T_e", ("core_boundary", "T_e")),
    ("core boundary T_i", ("core_boundary", "T_i")),
    ("core boundary pressure", ("core_boundary", "pressure_Pa")),
    ("electron content", ("electron_content",)),
    ("neutral content", ("neutral_content",)),
    ("stored energy", ("stored_energy_J",)),
    ("sheath power, electrons", ("power_W", "sheath_electrons")),
    ("sheath power, ions", ("power_W", "sheath_ions")),
    ("impurity radiation", ("power_W", "impurity_radiation", "*")),
    ("neutral energy to walls", ("power_W", "neutrals_to_walls")),
    ("ion flow to walls", ("particles_per_s", "ions_to_walls")),
    ("lower outer target power", ("targets", "lower_outer", "power_W")),
    ("lower inner target power", ("targets", "lower_inner", "power_W")),
    ("main chamber wall power", ("targets", "main_chamber", "power_W")),
]


def lookup(d, path):
    for k in path:
        if k == "*":
            return sum(d.values())
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def compare(out, reference):
    """Table of this run against a reference summary; returns the largest relative difference"""
    worst = 0.0
    print("%-34s %12s %12s %9s" % ("quantity", "reference", "this run", "rel. diff"))
    for label, path in COMPARE:
        a, b = lookup(reference, path), lookup(out, path)
        if a is None or b is None:
            continue
        rel = abs(b - a) / max(abs(a), 1e-300)
        worst = max(worst, rel)
        print("%-34s %12.4g %12.4g %8.1f%%" % (label, a, b, 100 * rel))
    return worst


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("run")
    p.add_argument("grid")
    p.add_argument("--json")
    p.add_argument("--figure")
    p.add_argument("--compare", help="reference summary.json to compare with")
    p.add_argument("--index", type=int, default=-1)
    p.add_argument("--gamma-e", type=float, default=4.5)
    p.add_argument("--gamma-i", type=float, default=2.5)
    a = p.parse_args()
    out, g, fields = summarise(a.run, a.grid, a.index, a.gamma_e, a.gamma_i)
    text = json.dumps(out, indent=1)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")
    brief = {k: v for k, v in out.items() if k not in ("outer_midplane", "history")}
    brief["outer_midplane_separatrix"] = out["outer_midplane"]["separatrix"]
    print(json.dumps(brief, indent=1))
    if a.compare:
        with open(a.compare) as f:
            worst = compare(out, json.load(f))
        print("largest relative difference: %.1f%%" % (100 * worst))
    if a.figure:
        figure(a.figure, g, fields, "%s, t = %.1f ms" % (a.run, 1e3 * out["time_s"]))


if __name__ == "__main__":
    main()
