"""Agreement tests on the TCV-X21 grid (run from us_runs/ after the timing runs).

usage: python numerics/agreement.py setup|analyse [v17|v18]
       (v18, the default: current bound on the skew correction; v17: the first bound,
        output in numerics/agreement)

One Euler step of dn/dt = div(D grad_perp n) alone (D = 1 m^2/s, skew
correction with its bound, as in the runs) from the TCV-X21 ion density, on
1, 2, 4 and 8 processes. dn/dt from the code is compared between the
decompositions and with tools/unstructured/usops.py, written independently
in numpy (serial).
(An explicit step of the full model is not possible from the TCV-X21 steady
state: the atomic rates are invalid in a few near-vacuum cells at the inner
wall, which the implicit solver never evaluates on their own.)
   The step is long enough (dt |dn/dt| comparable to n) that dn/dt is
   recovered from the two states to round-off.
"""
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, "../tools/unstructured")

VERSION = sys.argv[2] if len(sys.argv) > 2 else "v18"
BIN = "bin/%s/hermes-3" % VERSION
GRID = "grids/tcvx21_goat_drifts.grd.nc"
OUT = Path("numerics/agreement" if VERSION == "v17" else "numerics/agreement_" + VERSION)
# the bound on the skew correction that the binary uses by default (usops arguments)
BOUND = dict(limit=1.0, smooth=0.05, harmonic=False) if VERSION == "v17" else {}
OMEGA = 9.5788e7          # Omega_ci of the TCV-X21 runs [1/s]; read back from the dumps below
DIFFUSION_INPUT = """nout = 1
timestep = {dt}
MXG = 0
MZ = 1

[mesh]
file = "{grid}"

[solver]
type = euler

[hermes]
components = i
Nnorm = 1e19
Bnorm = 1
Tnorm = 50

[i]
type = evolve_density, anomalous_diffusion
AA = 2
charge = 1
anomalous_D = 1
core_density = 1.2e19

[Ni]
function = 1
bndry_all = neumann
"""


def run(args, log):
    with open(log, "w") as f:
        rc = subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, env={**__import__("os").environ,
                            "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}).returncode
    print(" ".join(args[:6]), "... exit", rc, flush=True)
    return rc


def setup():
    OUT.mkdir(parents=True, exist_ok=True)
    # 2. diffusion alone from the ion density of the steady state, on 1, 2, 4 and 8 processes
    for n in (1, 2, 4, 8):
        d = OUT / ("diffusion_np%d" % n)
        d.mkdir(exist_ok=True)
        for f in d.glob("BOUT.*"):
            f.unlink()
        (d / "BOUT.inp").write_text(DIFFUSION_INPUT.format(dt=1.0, grid=GRID))
        run(["mpirun", "-np", str(n), BIN, "-d", str(d), "nout=0"], d / "run0.log")   # restart files with the layout
        for k in range(n):
            src = netCDF4.Dataset("t22_timing/restart_np%d/BOUT.restart.%d.nc" % (n, k))
            with netCDF4.Dataset(str(d / ("BOUT.restart.%d.nc" % k)), "a") as dst:
                dst["Ni"][:] = src["Nd+"][:]
            src.close()
        for f in d.glob("BOUT.dmp.*"):
            f.unlink()
        run(["mpirun", "-np", str(n), BIN, "-d", str(d), "restart"], d / "run.log")


def analyse():
    import flatgrid
    import usdata
    import usops
    print("cross-field diffusion alone (skew correction with its bound), one Euler step on the TCV-X21 grid")
    rates = {}
    for n in (1, 2, 4, 8):
        d = str(OUT / ("diffusion_np%d" % n))
        r = usdata.Run(d, GRID)
        N = np.asarray(r("Ni"))
        with netCDF4.Dataset(d + "/BOUT.dmp.0.nc") as ds:
            t = np.asarray(ds["t_array"][:]) / float(ds["Omega_ci"][...])
        rates[n] = (N[1][:r.n] - N[0][:r.n]) / (t[1] - t[0])
    for n in (2, 4, 8):
        print("   %d processes against 1: largest difference %.1e of the largest rate (bit for bit: %s)"
              % (n, np.abs(rates[n] - rates[1]).max() / np.abs(rates[1]).max(), bool((rates[n] == rates[1]).all())))
    d = str(OUT / "diffusion_np8")
    r = usdata.Run(d, GRID)
    N = np.asarray(r("Ni"))
    with netCDF4.Dataset(d + "/BOUT.dmp.0.nc") as ds:
        t = np.asarray(ds["t_array"][:]) / float(ds["Omega_ci"][...])
    g = flatgrid.read(GRID)
    nreal = g["nreal"]
    ntot = len(g["vol"])
    f = np.zeros(ntot)
    f[:nreal] = N[0][:nreal]
    core = g["face_kind"] == flatgrid.KIND_CORE
    wall = g["face_kind"] == flatgrid.KIND_WALL
    f[g["face_R"][core]] = 1.2e19
    f[g["face_R"][wall]] = f[g["face_L"][wall]]
    ref_rate = usops.div_a_grad_perp_upwind(g, np.ones(ntot), f, isotropic=True, skew=True, **BOUND)[:nreal]
    code_rate = (N[1][:nreal] - N[0][:nreal]) / (t[1] - t[0])
    err = np.abs(code_rate - ref_rate)
    print("   dt = %.3e s; largest |dn/dt| %.3e m^-3/s; largest |code - numpy| %.2e of it; RMS %.2e"
          % (t[1] - t[0], np.abs(ref_rate).max(), err.max() / np.abs(ref_rate).max(),
             np.sqrt(np.mean(err**2)) / np.abs(ref_rate).max()))
    # the same without the correction, to show that the comparison is sensitive to it
    two = usops.div_a_grad_perp_upwind(g, np.ones(ntot), f, isotropic=True, skew=False)[:nreal]
    print("   (two-point reference instead: largest difference %.2e of the largest rate)"
          % (np.abs(code_rate - two).max() / np.abs(ref_rate).max()))
    # conservation: interior flows cancel and wall faces carry nothing, so the
    # volume integral of dn/dt is the flow in through the core faces
    grad = usops.face_gradient(g, f, True, skew=True, **BOUND)
    inflow = grad[core].sum()                    # coefficient D = 1 on every face
    total = (code_rate * g["vol"][:nreal]).sum()
    scale = (np.abs(code_rate) * g["vol"][:nreal]).sum()
    print("   conservation: volume integral of dn/dt %.6e /s/rad, flow in through the core %.6e; "
          "difference %.1e of the integral of |dn/dt|" % (total, inflow, abs(total - inflow) / scale))


if __name__ == "__main__":
    {"setup": setup, "analyse": analyse}[sys.argv[1]]()
