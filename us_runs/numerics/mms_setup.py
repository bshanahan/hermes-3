"""Method of manufactured solutions for the flat-cell code (Hermes-3 itself).

usage: python mms_setup.py <out dir> [resolutions, default 4,8,16,32,64]

Problem: steady cross-field diffusion of a density, dn/dt = div(D grad_perp n) + S,
with D = 1 m^2/s, the whole perpendicular gradient (Hermes-3's default on flat
cells), n fixed on the inner (core) boundary and no flow through the outer
wall. S is the exact cell average of div(-D grad_perp n_exact), from the
divergence theorem with Gauss quadrature on every edge, written into the grid
file as Ni_src. Grids: the skewed annulus of quadrilaterals or of triangles.

Manufactured solutions (s = (r - r_in) / w, both satisfy n = n_c at r_in and
dn/dr = 0 at r_out for every theta):
  smooth  n = n_c [1 + s(2 - s)(1 + 0.3 cos th + 0.2 sin 2th)]
  steep   n = n_c f(r) [1 + 0.3 s(2 - s) cos th],
          f = exp(-(r - r_in)/lam) + (r - r_in) exp(-w/lam)/lam, lam = 3 cm
          (falls by a factor 120 across the annulus, like an edge profile)
"""
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annuli  # noqa: E402
import flatgrid  # noqa: E402

N_C, D, LAM = 1e19, 1.0, 0.03
W = annuli.R_OUT - annuli.R_IN


def exact(name, x):
    """n and grad n (R, Z) at points x"""
    rad = x - np.array([annuli.R0, 0.0])
    r = np.hypot(*rad.T)
    th = np.arctan2(rad[:, 1], rad[:, 0])
    s = (r - annuli.R_IN) / W
    q, dq = s * (2 - s), (2 - 2 * s) / W                     # q(s) and dq/dr
    if name == "smooth":
        P = 1 + 0.3 * np.cos(th) + 0.2 * np.sin(2 * th)
        dP = -0.3 * np.sin(th) + 0.4 * np.cos(2 * th)
        n = N_C * (1 + q * P)
        dn_dr = N_C * dq * P
        dn_dth = N_C * q * dP
    elif name == "steep":
        e, eo = np.exp(-(r - annuli.R_IN) / LAM), np.exp(-W / LAM)
        f = e + (r - annuli.R_IN) * eo / LAM
        df = (-e + eo) / LAM
        M = 1 + 0.3 * q * np.cos(th)
        n = N_C * f * M
        dn_dr = N_C * (df * M + f * 0.3 * dq * np.cos(th))
        dn_dth = N_C * f * (-0.3 * q * np.sin(th))
    else:
        raise ValueError(name)
    e_r = rad / r[:, None]
    e_th = np.stack([-e_r[:, 1], e_r[:, 0]], 1)
    return n, dn_dr[:, None] * e_r + (dn_dth / r)[:, None] * e_th


def flux_density(name, x):
    """Gamma = -D grad_perp n (in-plane part) at points x"""
    _, gn = exact(name, x)
    b, _ = annuli.exact_unit_vectors(x)
    return -D * (gn - b * np.einsum("ij,ij->i", b, gn)[:, None])


def cell_source(g, name):
    """Cell average of S = div Gamma: outward flow of the exact Gamma through the faces / volume"""
    xg, wg = np.polynomial.legendre.leggauss(10)
    P = np.stack([g["vertex_R"], g["vertex_Z"]], 1)
    v1, v2 = g["face_v1"], g["face_v2"]
    n = np.stack([g["face_nR"], g["face_nZ"]], 1)          # from L to R
    flow = np.zeros(len(v1))
    for s, w in zip(0.5 * (xg + 1), 0.5 * wg):
        x = P[v1] + s * (P[v2] - P[v1])
        flow += w * np.einsum("ij,ij->i", flux_density(name, x), n) * x[:, 0]
    flow *= g["face_length"]                               # per radian, from L to R
    out = np.zeros(len(g["vol"]))
    np.add.at(out, g["face_L"], flow)
    np.add.at(out, g["face_R"], -flow)
    return out[:g["nreal"]] / g["vol"][:g["nreal"]]


INPUT = """# MMS: steady cross-field diffusion on a flat-cell annulus (written by mms_setup.py)
nout = 10
timestep = 95788 * 200     # 0.2 s; 2 s in all = 50 diffusion times
MXG = 0
MZ = 1

[mesh]
file = "{grid}"

[solver]
type = snes
pc_type = lu
max_nonlinear_iterations = 10
atol = 1e-12
rtol = 1e-10
stol = 1e-14
maxl = 260
matrix_free_operator = true
lag_jacobian = 50
timestep = 0.1
max_snes_failures = 30
scale_vars = true
timestep_control = pid_nonlinear_its
target_its = 4
kP = 0.7
kI = 0.3
kD = 0.0

[petsc]
snes_type = newtonls
ksp_type = fgmres
pc_factor_mat_solver_type = mumps
snes_linesearch_type = basic
snes_ksp_ew

[hermes]
components = i
Nnorm = 1e19
Bnorm = 1
Tnorm = 100

[i]
type = evolve_density, anomalous_diffusion
AA = 2
charge = 1
anomalous_D = 1            # m^2/s
core_density = 1e19        # m^-3, held at the core boundary

[Ni]
function = 1
bndry_all = neumann
"""


if __name__ == "__main__":
    out = Path(sys.argv[1])
    res = [int(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else [4, 8, 16, 32, 64]
    (out / "grids").mkdir(parents=True, exist_ok=True)
    for tri in (False, True):
        for n in res:
            for name in ("smooth", "steep"):
                tag = "%s_%s_n%d" % ("tri" if tri else "quad", name, n)
                path = out / "grids" / (tag + ".grd.nc")
                g = annuli.build(n, 4 * n, skew=0.4, triangles=tri)
                flatgrid.write(g, str(path), blocks=48)
                gr = flatgrid.read(str(path))
                src = cell_source(gr, name)
                s_of_l = gr["storage_index"]
                with netCDF4.Dataset(str(path), "a") as d:
                    a = np.zeros(len(d.dimensions["y"]))
                    a[s_of_l[:gr["nreal"]]] = src
                    v = d.createVariable("Ni_src", "f8", ("x", "y"))
                    v[0, :] = a
                    v.units = "m^-3 s^-1"
                    d.setncattr("mms", "mms_setup.py %s" % name)
                (out / (tag + ".inp")).write_text(INPUT.format(grid="grids/" + tag + ".grd.nc"))
                print(tag, gr["nreal"], "cells; source range %.3e .. %.3e" % (src.min(), src.max()), flush=True)
