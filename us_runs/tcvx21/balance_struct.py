"""Particle and power balance of a structured-grid TCV-X21 run (NXPE = 1).

usage (from us_runs/): python tcvx21/balance_struct.py <run dir>
"""
import sys, glob
import numpy as np, netCDF4

run = sys.argv[1]
files = sorted(glob.glob(run + "/BOUT.dmp.*.nc"), key=lambda f: int(f.split(".")[-2]))
names = ("Sd+_sheath", "Sd_pump", "Sd_target_recycle", "Sd_wall_recycle", "Ee_sheath", "Ed+_sheath", "Pd+_src", "Pe_src",
         "Nd+", "Nd", "Ed_pump", "Sd+_iz", "Rd+_ex", "Rd+_rec")
tot = {}
for f in files:
    with netCDF4.Dataset(f) as d:
        myg, mxg = int(d["MYG"][...]), int(d["MXG"][...])
        J, dx, dy = [np.asarray(d[k][:]) for k in ("J", "dx", "dy")]
        J, dx, dy = [(a[..., 0] if a.ndim == 3 else a) for a in (J, dx, dy)]
        dV = (2*np.pi*J*dx*dy)[mxg:-mxg, myg:-myg]*float(d["rho_s0"][...])**3      # m^3 (J dx is normalised to rho_s0^3)
        t = np.asarray(d["t_array"][:])/float(d["Omega_ci"][...])
        for name in names:
            if name not in d.variables: continue
            v = d[name]; c = v.getncattr("conversion") if "conversion" in v.ncattrs() else 1.0
            for it, tag in ((-1, ""), (-6, " before")):
                tot[name + tag] = tot.get(name + tag, 0.0) + (np.asarray(v[it, mxg:-mxg, myg:-myg, 0])*c*dV).sum()
g = lambda k: tot.get(k, 0.0)
dt = t[-1] - t[-6]
print("t = %.1f ms | ion content %.3e, neutral content %.3e | relative change per 10 ms: ions %.1e, neutrals %.1e" % (
    1e3*t[-1], g("Nd+"), g("Nd"), (g("Nd+")/g("Nd+ before") - 1)*0.01/dt, (g("Nd")/g("Nd before") - 1)*0.01/dt))
print("particles [1/s]: ion flow to targets %.3e | recycled at targets %.3e | recycled at radial edges %.3e | neutrals absorbed at radial edges %.3e | ionisation %.3e" % (
    -g("Sd+_sheath"), g("Sd_target_recycle"), g("Sd_wall_recycle"), -g("Sd_pump"), g("Sd+_iz")))
P_in = 1.5*(g("Pd+_src") + g("Pe_src")); P_sh = -(g("Ee_sheath") + g("Ed+_sheath"))
print("power [kW]: heating %.1f | sheath %.1f (electrons %.1f, ions %.1f) | 13.6 eV per ion %.1f | radiation %.1f | remainder %.1f" % (
    P_in/1e3, P_sh/1e3, -g("Ee_sheath")/1e3, -g("Ed+_sheath")/1e3, 13.6*1.602e-19*(-g("Sd+_sheath"))/1e3, -(g("Rd+_ex") + g("Rd+_rec"))/1e3, (P_in - P_sh)/1e3))
