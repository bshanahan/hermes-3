"""Orthogonal (hypnotoad) grid for the TCV-X21 equilibrium, for comparison with the GOAT grid.

usage: python -I make_grid.py <TCV-X21 eqdsk> <out dir> [key=value ...]
   e.g. python -I make_grid.py ../data/1.experimental_data_reference_scenario_65402_t1.eqdsk . nx_sol=22

Two things are changed in the equilibrium file before hypnotoad reads it:
  * the poloidal flux is divided by 2 pi. The TCV-X21 file stores psi in Wb, hypnotoad
    expects Wb/rad. Checked with Ampere's law against the plasma current in the header.
  * the limiter is replaced by the 8-segment TCV wall that the GOAT grid uses, so that
    both grids end on the same surfaces.
"""
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from hypnotoad.geqdsk._geqdsk import read as read_geqdsk, write as write_geqdsk
from hypnotoad.cases import tokamak
from hypnotoad.core.mesh import BoutMesh

src, out = sys.argv[1], sys.argv[2].rstrip("/")
overrides = dict(a.split("=", 1) for a in sys.argv[3:])

with open(src) as f:
    g = read_geqdsk(f)

# ---- flux units (Ampere's law on the box boundary)
R = g["rleft"] + g["rdim"]*np.arange(g["nx"])/(g["nx"] - 1)
Z = g["zmid"] - g["zdim"]/2 + g["zdim"]*np.arange(g["ny"])/(g["ny"] - 1)
spl = RectBivariateSpline(R, Z, np.asarray(g["psi"]))
m = 4000
def side(r0, z0, r1, z1):
    r = np.linspace(r0, r1, m); z = np.linspace(z0, z1, m)
    BR = -spl(r, z, dy=1, grid=False)/r; BZ = spl(r, z, dx=1, grid=False)/r
    return np.trapz(BR, r) + np.trapz(BZ, z)
e = 2e-3
circ = (side(R[0]+e, Z[0]+e, R[-1]-e, Z[0]+e) + side(R[-1]-e, Z[0]+e, R[-1]-e, Z[-1]-e)
        + side(R[-1]-e, Z[-1]-e, R[0]+e, Z[-1]-e) + side(R[0]+e, Z[-1]-e, R[0]+e, Z[0]+e))
ratio = abs(circ/(4e-7*np.pi)/g["cpasma"])
if abs(ratio - 1) < 0.05:
    div = 1.0
elif abs(ratio/(2*np.pi) - 1) < 0.05:
    div = 2*np.pi
else:
    raise SystemExit("cannot tell the flux units: enclosed current / header current = %.3f" % ratio)
print("enclosed current / header current = %.3f -> psi divided by %.4f" % (ratio, div))
for k in ("psi", "simagx", "sibdry"):
    g[k] = np.asarray(g[k])/div if k == "psi" else g[k]/div
for k in ("pprime", "ffprime"):
    g[k] = np.asarray(g[k])*div

# ---- the wall of the GOAT grid (closed, anticlockwise in the R-Z plane)
wall = np.array([(0.671, 0.750), (0.624, 0.700), (0.624, -0.700), (0.671, -0.750),
                 (0.965, -0.750), (1.136, -0.543), (1.136, 0.543), (0.965, 0.750)])
g["rlim"], g["zlim"] = wall[:, 0], wall[:, 1]
g["nlim"] = len(wall)
eq_path = out + "/tcvx21_psi_per_rad.geqdsk"
with open(eq_path, "w") as f:
    write_geqdsk(g, f)

settings = dict(
    orthogonal=True, shiftedmetric=True, y_boundary_guards=0, number_of_processors=1,
    curvature_type="curl(b/B)",
    psinorm_core=0.90,            # as the GOAT grid
    psinorm_sol=1.24,             # the inner wall touches psiN = 1.257 above the X-point
    psinorm_pf=0.75,
    nx_core=10, nx_sol=22, nx_pf=10,
    ny_inner_divertor=16, ny_sol=64, ny_outer_divertor=24,
    finecontour_Nfine=100, refine_width=1e-5, refine_atol=2e-8,
    xpoint_poloidal_spacing_length=0.05, target_all_poloidal_spacing_length=0.05,
)
for k, v in overrides.items():
    settings[k] = type(settings.get(k, 0.0))(v) if not isinstance(settings.get(k), bool) else v == "True"
name = "%s/tcvx21_%dx%d.grd.nc" % (out, settings["nx_core"] + settings["nx_sol"],
                                    settings["ny_inner_divertor"] + settings["ny_sol"] + settings["ny_outer_divertor"])
with open(eq_path) as f:
    eq = tokamak.read_geqdsk(f, settings=settings, nonorthogonal_settings=settings)
mesh = BoutMesh(eq, settings)
mesh.geometry()
mesh.writeGridfile(name)
print("wrote", name)
