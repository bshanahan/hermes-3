"""Write GOAT's equilibrium input (rzpsi.dat) from a G-EQDSK file.

usage: python geqdsk_to_rzpsi.py equilibrium.geqdsk rzpsi.dat
           [--box RMIN RMAX ZMIN ZMAX] [--points NR NZ]
           [--psi-units auto|wb_per_rad|wb]

GOAT reads the poloidal flux psi in Wb on a uniform (R, Z) grid. This script

  * resamples psi with a bicubic spline onto the requested box and number of
    points (default: the box and points of the G-EQDSK file). The box must
    contain the whole vessel given to GOAT and should leave out the coils;
  * converts psi to Wb. A G-EQDSK file normally holds psi in Wb/rad, but
    some hold Wb. The units are found with Ampere's law: the current enclosed
    by the plasma boundary, computed from the psi map as if it were in Wb/rad,
    must equal the plasma current in the header (ratio 1), or be 2 pi times
    larger (file in Wb). Anything else is an error;
  * shifts psi so that the separatrix is at zero and gives it the sign that
    makes it increase outwards from the magnetic axis;
  * writes psi with 10 significant digits. With 6 digits GOAT's topology
    construction fails near the X-points ("lines not allocated for cell").

The last lines printed are

    RBTOR=<R * B_toroidal at the plasma edge, T m>

which is the value for 'goat.mf.RBtor' in GOAToptions.dat, and the box, to be
compared with the vessel in GOAT's structure.dat.

Needs numpy and scipy only.
"""
import argparse
import re
import numpy as np
from scipy.interpolate import RectBivariateSpline

NUMBER = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?")


def read_geqdsk(path):
    """Minimal G-EQDSK reader: header, fpol, psi(R, Z) and the plasma boundary."""
    with open(path) as f:
        header = f.readline()
        body = f.read()
    nw, nh = (int(v) for v in header.split()[-2:])
    numbers = NUMBER.findall(body)
    pos = 0

    def take(n):
        nonlocal pos
        out = np.array([float(v.replace("D", "E").replace("d", "e")) for v in numbers[pos:pos + n]])
        if out.size != n:
            raise ValueError("%s ends early: wanted %d numbers at position %d" % (path, n, pos))
        pos += n
        return out

    rdim, zdim, rcentr, rleft, zmid = take(5)
    rmaxis, zmaxis, simagx, sibdry, bcentr = take(5)
    cpasma = take(5)[0]
    take(5)
    fpol = take(nw)
    take(3 * nw)                            # pressure, FF', p'
    psi = take(nw * nh).reshape(nh, nw).T   # file has R varying fastest; psi[R index, Z index]
    take(nw)                                # q
    rbdry = zbdry = np.zeros(0)
    if pos + 2 <= len(numbers):
        nbdry, nlim = (int(float(v)) for v in numbers[pos:pos + 2])
        pos += 2
        boundary = take(2 * nbdry)
        rbdry, zbdry = boundary[0::2], boundary[1::2]
    return dict(nw=nw, nh=nh, rdim=rdim, zdim=zdim, rcentr=rcentr, rleft=rleft, zmid=zmid,
                rmaxis=rmaxis, zmaxis=zmaxis, simagx=simagx, sibdry=sibdry, bcentr=bcentr,
                cpasma=cpasma, fpol=fpol, psi=psi, rbdry=rbdry, zbdry=zbdry)


def psi_units(g, spline):
    """Factor that turns the file's psi into Wb/rad (1 or 1/(2 pi)), from Ampere's law."""
    Rb, Zb = g["rbdry"], g["zbdry"]
    if Rb.size < 8 or g["cpasma"] == 0:
        raise SystemExit("The file has no plasma boundary or no plasma current, so the units of psi "
                         "cannot be checked. Pass --psi-units wb_per_rad (usual) or --psi-units wb.")
    Rm, Zm = 0.5 * (Rb + np.roll(Rb, -1)), 0.5 * (Zb + np.roll(Zb, -1))
    dR, dZ = np.roll(Rb, -1) - Rb, np.roll(Zb, -1) - Zb
    BR = -spline(Rm, Zm, dy=1, grid=False) / Rm
    BZ = spline(Rm, Zm, dx=1, grid=False) / Rm
    ratio = abs(np.sum(BR * dR + BZ * dZ) / (4e-7 * np.pi) / g["cpasma"])
    if abs(ratio - 1) < 0.05:
        print("psi in the file is in Wb/rad (enclosed current / header current = %.3f)" % ratio)
        return 1.0
    if abs(ratio / (2 * np.pi) - 1) < 0.05:
        print("psi in the file is in Wb, not Wb/rad (enclosed current / header current = %.3f = 2 pi)" % ratio)
        return 1.0 / (2 * np.pi)
    raise SystemExit("psi map and plasma current disagree: enclosed / header current = %.3f "
                     "(expected 1 for Wb/rad or 6.283 for Wb)" % ratio)


def block(values, fmt):
    values = np.asarray(values).ravel()
    rows = [values[i:i + 6] for i in range(0, values.size, 6)]
    return "\n".join(" " + "".join(fmt % v for v in row) for row in rows) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("geqdsk")
    p.add_argument("rzpsi")
    p.add_argument("--box", type=float, nargs=4, metavar=("RMIN", "RMAX", "ZMIN", "ZMAX"))
    p.add_argument("--points", type=int, nargs=2, metavar=("NR", "NZ"))
    p.add_argument("--psi-units", choices=("auto", "wb_per_rad", "wb"), default="auto")
    a = p.parse_args()

    g = read_geqdsk(a.geqdsk)
    r = g["rleft"] + g["rdim"] * np.arange(g["nw"]) / (g["nw"] - 1)
    z = g["zmid"] - g["zdim"] / 2 + g["zdim"] * np.arange(g["nh"]) / (g["nh"] - 1)
    spline = RectBivariateSpline(r, z, g["psi"], kx=3, ky=3)

    per_rad = {"auto": None, "wb_per_rad": 1.0, "wb": 1.0 / (2 * np.pi)}[a.psi_units]
    if per_rad is None:
        per_rad = psi_units(g, spline)

    box = a.box or (r[0], r[-1], z[0], z[-1])
    if box[0] < r[0] - 1e-9 or box[1] > r[-1] + 1e-9 or box[2] < z[0] - 1e-9 or box[3] > z[-1] + 1e-9:
        raise SystemExit("--box %s lies outside the equilibrium (R %.4f..%.4f, Z %.4f..%.4f)"
                         % (box, r[0], r[-1], z[0], z[-1]))
    nr, nz = a.points or (g["nw"], g["nh"])
    R = np.linspace(box[0], box[1], nr)
    Z = np.linspace(box[2], box[3], nz)

    sign = -1.0 if g["simagx"] > g["sibdry"] else 1.0
    to_weber = 2 * np.pi * per_rad
    psi = sign * to_weber * (spline(R, Z) - g["sibdry"])   # [nr, nz], Wb, separatrix at 0

    with open(a.rzpsi, "w") as f:
        f.write(" \n %d  %d\n $r\n nr=%d\n" % (nr, nz, nr))
        f.write(block(R, "%15.7E"))
        f.write(" \n $z\n nz=%d\n" % nz)
        f.write(block(Z, "%15.7E"))
        f.write(" \n $psi\n")
        f.write(block(psi.T, "%18.9E"))   # R varies fastest

    print("wrote %s: %d x %d points" % (a.rzpsi, nr, nz))
    print("psi on the magnetic axis %.6f Wb, separatrix 0 Wb (sign applied %+d)"
          % (sign * to_weber * (g["simagx"] - g["sibdry"]), sign))
    print("BOX=%.4f %.4f %.4f %.4f" % tuple(box))
    print("RBTOR=%.10g" % g["fpol"][-1])


if __name__ == "__main__":
    main()
