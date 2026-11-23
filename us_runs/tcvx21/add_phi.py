"""Add a first guess of the potential to restart files: phi = 3 T_e in every cell.

usage: python tcvx21/add_phi.py <directory with BOUT.restart.*.nc>
The sheath potential of a floating deuterium plasma is close to 3 T_e.
"""
import glob
import sys
import netCDF4

for path in sorted(glob.glob(sys.argv[1] + "/BOUT.restart.*.nc")):
    with netCDF4.Dataset(path, "a") as d:
        ne = d["Nd+"][:]
        te = d["Pe"][:] / ne
        v = d["phi"] if "phi" in d.variables else d.createVariable("phi", d["Pe"].dtype, d["Pe"].dimensions)
        for a in d["Pe"].ncattrs():
            v.setncattr(a, d["Pe"].getncattr(a))
        v[:] = 3.0 * te
    print(path, "phi added")
