"""Write GOAT's rzpsi.dat from a geqdsk file, on the geqdsk's own box.

usage: python make_rzpsi.py equilibrium.geqdsk rzpsi.dat [nr nz]

GOAT's Fortran reader takes psi in Wb, so a geqdsk value in Wb/rad is
multiplied by 2*pi. psi is shifted so the separatrix is at zero, given the
sign that makes it increase outwards from the magnetic axis (as in the SPARC
case), and written with 10 significant digits.

The units of psi in the file are checked with Ampere's law: the current
enclosed by the plasma boundary, computed from the psi map as if it were in
Wb/rad, must equal the plasma current in the header. Some files store psi in
Wb instead (the TCV-X21 reference equilibrium does); the computed current is
then 2*pi too large, and psi is used as it is. Anything else is an error.
"""
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from hypnotoad.geqdsk._geqdsk import read

geq, out = sys.argv[1:3]
with open(geq) as f:
    g = read(f)
r = g["rleft"] + g["rdim"] * np.arange(g["nx"]) / (g["nx"] - 1)
z = g["zmid"] - g["zdim"] / 2 + g["zdim"] * np.arange(g["ny"]) / (g["ny"] - 1)
spline = RectBivariateSpline(r, z, np.asarray(g["psi"]), kx=3, ky=3)
nr, nz = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (g["nx"], g["ny"])
R = np.linspace(r[0], r[-1], nr)
Z = np.linspace(z[0], z[-1], nz)
# Units of psi from Ampere's law around the plasma boundary
Rb, Zb = np.asarray(g["rbdry"]), np.asarray(g["zbdry"])
Rm, Zm = 0.5 * (Rb + np.roll(Rb, -1)), 0.5 * (Zb + np.roll(Zb, -1))
dR, dZ = np.roll(Rb, -1) - Rb, np.roll(Zb, -1) - Zb
BR = -spline(Rm, Zm, dy=1, grid=False) / Rm
BZ = spline(Rm, Zm, dx=1, grid=False) / Rm
ratio = abs(np.sum(BR * dR + BZ * dZ) / (4e-7 * np.pi) / g["cpasma"])
if abs(ratio - 1) < 0.05:
    to_weber = 2 * np.pi
    print("psi in the file is in Wb/rad (enclosed current / header current = %.3f)" % ratio)
elif abs(ratio / (2 * np.pi) - 1) < 0.05:
    to_weber = 1.0
    print("psi in the file is in Wb, not Wb/rad (enclosed current / header current = %.3f = 2 pi)" % ratio)
else:
    raise SystemExit("psi map and plasma current disagree: enclosed / header current = %.3f" % ratio)

sign = -1.0 if g["simagx"] > g["sibdry"] else 1.0
psi = sign * to_weber * (spline(R, Z) - g["sibdry"])  # [nr, nz], Wb, separatrix at 0


def block(a, fmt="%13.5E"):
    a = np.asarray(a).ravel()
    rows = [a[i:i + 6] for i in range(0, a.size, 6)]
    return "\n".join(" " + "".join(fmt % v for v in row) for row in rows) + "\n"


with open(out, "w") as f:
    f.write(" \n %d  %d\n $r\n nr=%d\n" % (R.size, Z.size, R.size))
    f.write(block(R, "%15.7E"))
    f.write(" \n $z\n nz=%d\n" % Z.size)
    f.write(block(Z, "%15.7E"))
    f.write(" \n $psi\n")
    f.write(block(psi.T, "%18.9E"))  # r varies fastest
print("R*Btor at the edge [T m]: %.10g" % g["fpol"][-1])
print("psi axis %.6f Wb, separatrix 0 Wb, sign applied %+d" % (sign * to_weber * (g["simagx"] - g["sibdry"]), sign))
