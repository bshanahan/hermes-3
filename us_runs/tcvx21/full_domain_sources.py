"""Heating and fuelling sources for a TCV-X21 run on a GOAT grid without a core boundary.

usage (from us_runs/):
    python tcvx21/full_domain_sources.py <grid in> <grid out> <particles per second>
        [power W, default 160e3] [heated inside psiN, default 0.8]
        [fuelled between psiN, default "0.8,0.9"]

Writes a copy of the grid file with three per-cell arrays that Hermes-3 reads
as sources (evolve_pressure / evolve_density read P<species>_src [Pa/s] and
N<species>_src [m^-3 s^-1] from the mesh):
    Pd+_src, Pe_src   half of the power each, uniform per volume in the closed
                      flux tubes inside the given psiN
    Nd+_src           the particle source, uniform per volume in the closed
                      flux tubes between the two psiN values
psiN is the flux label of each cell from the grid file (the mean psi of its
vertices), so a flux tube is either wholly in a region or wholly out of it.
The source particles carry no energy: all of the power is in the heating.
"""
import shutil
import sys

import netCDF4
import numpy as np

sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import flatgrid
import compare as c


def main(grid_in, grid_out, particles, power=160e3, heat_psin=0.8, fuel=(0.8, 0.9)):
    g = flatgrid.read(grid_in)
    n = g["nreal"]
    twopi = 2 * np.pi
    psin = 1 + g["psi"][:n] * twopi / c.q["simagx"]     # as compare.load: grid psi in Wb/rad, sign flipped
    tube = g["tube"][:n]
    closed_tube = np.array([(g["prev"][:n][tube == t] >= 0).all() and (g["next"][:n][tube == t] >= 0).all()
                            for t in range(tube.max() + 1)])
    closed = closed_tube[tube]
    vol = twopi * g["vol"][:n]                            # m^3
    heated = closed & (psin < heat_psin)
    fuelled = closed & (psin > fuel[0]) & (psin < fuel[1])
    if not heated.any() or not fuelled.any():
        raise SystemExit("empty source region: heated %d cells, fuelled %d cells" % (heated.sum(), fuelled.sum()))
    Vh, Vf = vol[heated].sum(), vol[fuelled].sum()
    p_src = np.where(heated, (2.0 / 3.0) * 0.5 * power / Vh, 0.0)   # Pa/s for each of ions and electrons
    n_src = np.where(fuelled, particles / Vf, 0.0)
    print("heating %.1f kW (half to ions, half to electrons) in %d cells, %.3f m^3, psiN < %.2f (flux tubes %s)" % (
        power / 1e3, heated.sum(), Vh, heat_psin, np.unique(tube[heated]).size))
    print("fuelling %.3e /s in %d cells, %.3f m^3, %.2f < psiN < %.2f" % (particles, fuelled.sum(), Vf, *fuel))
    print("check: power %.1f kW, particles %.3e /s" % (1.5 * 2 * (p_src * vol).sum() / 1e3, (n_src * vol).sum()))

    shutil.copyfile(grid_in, grid_out)
    s_of_l = g["storage_index"]
    with netCDF4.Dataset(grid_out, "a") as d:
        N = len(d.dimensions["y"])
        for name, values, units in (("Pd+_src", p_src, "Pa/s"), ("Pe_src", p_src, "Pa/s"), ("Nd+_src", n_src, "m^-3 s^-1")):
            a = np.zeros(N)
            a[s_of_l[:n]] = values
            v = d.createVariable(name, "f8", ("x", "y"))
            v[0, :] = a
            v.units = units
        d.setncattr("sources", "full_domain_sources.py: %.1f kW heating inside psiN %.2f, %.3e particles/s between psiN %.2f and %.2f"
                    % (power / 1e3, heat_psin, particles, *fuel))
    print("wrote", grid_out)


if __name__ == "__main__":
    args = sys.argv[1:]
    fuel = tuple(float(v) for v in args[5].split(",")) if len(args) > 5 else (0.8, 0.9)
    main(args[0], args[1], float(args[2]), float(args[3]) if len(args) > 3 else 160e3,
         float(args[4]) if len(args) > 4 else 0.8, fuel)
