# Hermes-3 on a wall-conforming SPARC grid: set-up and sample run

These instructions take a machine with nothing installed to a finished edge
transport run of a SPARC L-mode plasma, and say how to drive that run from an
integrated-modelling framework. They are written to be followed step by step by
an automated agent. Every command below was run from an empty directory on
Ubuntu 24.04 on 24 November 2026; the times quoted are from that run on 8 cores
of a workstation that was busy with other jobs.

What you get:

| Piece | What it is | Where it comes from |
| --- | --- | --- |
| Hermes-3 with flat-cell ("unstructured") grids | Plasma edge fluid code (BOUT++). This branch adds finite-volume operators on cells of any shape | `https://github.com/bshanahan/hermes-3`, branch `unstructured` |
| GOAT | Mesh generator: field-aligned cells that reach the real wall | `https://github.com/TME-TFE/goat`, branch `release/1.0.0-beta` |
| Equilibrium | SPARC 12.2 T, 8.5 MA L-mode, lower single null, "VV" strike points (G-EQDSK) | `https://github.com/cfs-energy/SPARCPublic`, `LmodeLowerSingleNull/VV` |

The chain is

    equilibrium.geqdsk --geqdsk_to_rzpsi.py--> rzpsi.dat --GOAT--> traduit.out.b2us
        --goat_to_bout.py--> sparc_lmode.grd.nc --Hermes-3--> BOUT.dmp.*.nc --summary.py--> summary.json

The branch is experimental research code. It is not part of upstream Hermes-3,
and the grid format, the options and the results may change.

## 1. Requirements

Verified with: Ubuntu 24.04, GCC and gfortran 13.3, CMake 3.28, Open MPI 4.1.6,
PETSc 3.19.6 built with MUMPS (the Ubuntu package), NetCDF 4.9.2 with the C++
interface 4.3.1, FFTW 3.3.10, Python 3.12.

    sudo apt install -y build-essential gfortran git cmake curl tcsh \
        libopenmpi-dev openmpi-bin petsc-dev libnetcdf-dev libnetcdf-c++4-dev libfftw3-dev \
        libopenblas-dev liblapack-dev libblas-dev libsuitesparse-dev python3-venv

This list is what the build uses on the verified machine; it was not replayed on a
freshly installed system.

On another system you need the same things by other means. Two points are not
optional:

- PETSc must include the MUMPS direct solver (the sample input uses
  `pc_factor_mat_solver_type = mumps`). A PETSc you build yourself needs
  `--download-mumps --download-scalapack`; then export `PETSC_DIR` and
  `PETSC_ARCH` before running CMake.
- GOAT needs `tcsh`, a Fortran 2008 compiler, OpenBLAS and SuiteSparse (UMFPACK, CXSparse).

Disk: 2 GB for the three source trees, the builds and the Python environment.
Memory: the sample run uses about 0.3 GB per process.

## 2. Download

Choose an empty working directory. All later paths are relative to it.

    export WORK=$PWD
    git clone --recursive -b unstructured https://github.com/bshanahan/hermes-3.git hermes-3-unstructured
    git clone -b release/1.0.0-beta https://github.com/TME-TFE/goat.git goat
    git -C goat checkout e0d3eefa291695f1f0130d4b1f5dff66d2b2aac5
    mkdir -p sparc && curl -L -o sparc/equilibrium.geqdsk \
        https://raw.githubusercontent.com/cfs-energy/SPARCPublic/5b913b8216d05346b22e2e1e87c2ff100761613f/LmodeLowerSingleNull/VV/equilibrium.geqdsk
    python3 -m venv venv && venv/bin/pip install numpy scipy netCDF4 matplotlib

Checks:

    git -C hermes-3-unstructured/external/BOUT-dev rev-parse HEAD
    # 6937d48a26bb5ddd3d9be1190f4a31576ee6ccd3
    sha256sum sparc/equilibrium.geqdsk
    # bf6aed3ff09488a96f566225a4e363c9237252f649ad8f370761e61eba4f846f

The BOUT++ commit lives on the branch `snes-jacobian-reuse` of
`boutproject/BOUT-dev`. If the recursive clone cannot find it, that branch has
been removed upstream: stop and report it, do not substitute another commit.

The Python tools need only numpy, scipy and netCDF4 (matplotlib for figures).

## 3. Build Hermes-3 (4 minutes on 8 cores)

    cd $WORK/hermes-3-unstructured
    git -C external/BOUT-dev diff --quiet && \
        git -C external/BOUT-dev apply ../../tools/unstructured/BOUT-dev.patch
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCHECK=0 \
        -DBOUT_USE_PETSC=ON -DBOUT_USE_SUNDIALS=OFF \
        -DHERMES_SLOPE_LIMITER=VanAlbada -DHERMES_UPDATE_GIT_SUBMODULE=OFF
    cmake --build build -j 8

- The patch changes two files of BOUT++ (the Jacobian pattern of flat-cell grids
  and a failure diagnostic). The first line applies it only if BOUT++ is unmodified,
  so the block can be repeated safely.
- `-DHERMES_UPDATE_GIT_SUBMODULE=OFF` is required. Without it the build resets
  BOUT++ to the recorded commit and silently discards the patch.
- `-DHERMES_SLOPE_LIMITER=VanAlbada` is the limiter all flat-cell results were
  obtained with (the default, MC, needs about twice as many time steps).
- CMake prints `ModuleNotFoundError: No module named 'cython'` and a warning
  about the BOUT++ Python interface. Ignore both.

Check the build (8 seconds; runs Hermes-3 on 1 and 3 processes against a
reference solution):

    HERMES_EXE=$PWD/build/hermes-3 $WORK/venv/bin/python tests/integrated/unstructured-diffusion/runtest
    # last line: " => Test passed"

The executable is `$WORK/hermes-3-unstructured/build/hermes-3`. It reads the
atomic rate tables from `$WORK/hermes-3-unstructured/json_database`, a path fixed
at compile time. If the source tree moves, add
`json_database_dir=/new/path/json_database` to every Hermes-3 command line.

## 4. Build GOAT (5 minutes, one core)

    cd $WORK/goat
    tcsh -c 'source setup.csh; make goat'
    ls executables/goat.exe

`setup.csh` must be sourced from the top directory of GOAT, as here. Warnings
about missing dependencies from `dpkg` and "Permission denied" lines from `find`
are harmless if `goat.exe` appears.

## 5. Make the grid (2 minutes)

    cd $WORK/sparc
    PYTHON=$WORK/venv/bin/python GOAT_THREADS=4 \
        $WORK/hermes-3-unstructured/examples/unstructured/make_grid.sh \
        $WORK/goat equilibrium.geqdsk sparc_lmode.grd.nc goat_case

The script (read its header) converts the equilibrium to GOAT's input, runs GOAT
with the vessel and meshing options of GOAT's own SPARC example, and converts
the mesh. Expected output, apart from paths:

    psi in the file is in Wb/rad (enclosed current / header current = 1.000)
    psi on the magnetic axis -14.990037 Wb, separatrix 0 Wb (sign applied +1)
    RBTOR=-22.496
    ...
    3402 cells + 314 ghost cells, 6817 faces
       ...
       parallel / cross-field transmissibility on flux-surface faces (max; must be < 1e-12) 0.0
       boundary faces core / wall                                         (75, 239)
    wrote .../sparc_lmode.grd.nc

GOAT is not reproducible to the cell: another machine, compiler or thread count
can give a slightly different mesh (3,407 cells in an earlier run here). A cell
count within a few of 3,402 and a script that ends with `wrote ...` mean success.
`goat_to_bout.py` stops with an assertion if the mesh is unusable.

The grid of the verified run is in the repository as
`examples/unstructured/sparc_lmode.grd.nc`. Use it to skip this step, or to
compare against your own grid.

What is in the grid: cells between the flux surface at normalised poloidal flux
0.90 (the "core boundary") and the wall, on both sides of the separatrix, with the
private-flux regions. Cells are stored as one long row (`nx = 1`), each followed
by ghost cells for its boundary faces. All geometry is in SI units per radian of
toroidal angle. `tools/unstructured/flatgrid.py` documents every variable.

## 6. The sample case

`examples/unstructured/sparc-lmode/BOUT.inp` is a complete input. Its header
names the inputs a framework would set:

| Option | Sample value | Meaning |
| --- | --- | --- |
| `power_core` | `17e6` | Heating power crossing the core boundary [W] |
| `ion_power_fraction` | `0.6` | Share of it carried by the ions (half D+, half T+); the rest goes to the electrons |
| `n_core` | `8e19` | Electron density held on the core boundary [m^-3], half D+ and half T+ |
| `neon_fraction` | `0.005` | Neon density over electron density, everywhere; sets the neon radiation |
| `D_perp` | `0.3` | Anomalous particle diffusivity [m^2/s], uniform |
| `chi_perp` | `0.15` | Anomalous heat diffusivity of ions and electrons [m^2/s], uniform |
| `recycle` | `1` | Fraction of the ions reaching any wall face that returns as neutrals |
| `pump_albedo` | `1` | Fraction of the neutrals reaching the pump surface that returns (1 = no pump) |
| `puff_rate` | `0` | Gas puff [atoms/s], half D and half T |

The rest of the model: D+ and T+ with their own densities, parallel momenta and
pressures; electrons with their own pressure, quasi-neutral, no net current; D and
T fluid neutrals (diffusive across the field, with parallel momentum) coupled by
ionisation, recombination and charge exchange; Braginskii parallel conduction
without a heat-flux limiter; a sheath at every wall face the field crosses; no
drifts. With `recycle = 1`, no pump and no puff, no particles cross the core
boundary in steady state, so `power_core` is all conducted.

Pump and valve: the input has a pump on the private-flux-region wall under the
lower X-point (the "dome", 1.7 m^2, `pump_boxes` in `[recycling]`) and a gas
valve at the outboard midplane (`locations` in `[puff_d]` and `[puff_t]`). Both
positions are generic placeholders, not SPARC's, and both are off by default.
With the core density held fixed, every pumped particle is replaced through the
core boundary and carries about 5/2 (T_e + T_i) of keV plasma with it. In a
detached divertor (neutral density about 3e21 m^-3 at the dome, seen on the way to
the stored state) a pump absorbing 1 % of the neutrals it meets
(`pump_albedo = 0.99`) would need 1e23 particles/s, far beyond what 17 MW can
carry. Keep `pump_albedo` at 0.9999 or closer to 1, or balance the pump with the
puff.

`core_power` is the total power through the core boundary: the code subtracts
what particle diffusion carries. The temperature on the boundary is a result.

## 7. Run

Always set the thread variables: without them OpenBLAS starts a thread per core
in every MPI process.

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    cd $WORK/hermes-3-unstructured/examples/unstructured

### 7a. The sample run (under 3 minutes on 8 processes)

`sparc-lmode/` holds the input and the converged state of the sample as restart
files for 8 processes on `sparc_lmode.grd.nc` of this directory. Copy the
directory (a run overwrites its restart files) and continue from it:

    cp -r sparc-lmode $WORK/runs/sample
    mpirun -np 8 $WORK/hermes-3-unstructured/build/hermes-3 -d $WORK/runs/sample restart \
        mesh:file=$PWD/sparc_lmode.grd.nc > $WORK/runs/sample/run.log 2>&1
    $WORK/venv/bin/python summary.py $WORK/runs/sample sparc_lmode.grd.nc \
        --json $WORK/runs/sample/summary.json --compare sparc-lmode/reference_summary.json

This runs 300 ms of plasma time (60 outputs, 5 ms apart). Measured: 2 min 38 s on
8 cores of a shared workstation. The last line of `summary.py` must read
`largest relative difference: 0.0%` (up to a few tenths of a percent on another
machine): the run stays on the stored steady state.

The input names its grid relative to the directory the command starts from;
`mesh:file` with a full path makes that irrelevant.

### 7b. Another grid, process count or physics: run_warm.sh

Restart files belong to one grid and one process count. `run_warm.sh` maps the
stored state onto any grid and process count and runs from it:

    PYTHON=$WORK/venv/bin/python NP=4 GRID=$WORK/sparc/sparc_lmode.grd.nc \
        ./run_warm.sh $WORK/runs/mine power_core=20e6 nout=40 timestep=957880

Options after the run directory go to Hermes-3 (`nout` outputs, `timestep` apart
in units of 1/95788 ms; here 40 outputs 10 ms apart). It ends with the main
numbers and writes `<run directory>/summary.json`. What it does, so that you can
do it by hand:

1. copies `sparc-lmode/BOUT.inp` into the run directory and runs Hermes-3 for
   one tiny step, only to create restart files for this grid and process count;
2. overwrites the evolving fields in them with the stored state, using
   `tools/unstructured/regrid_restart.py`: each cell takes the value of a nearby
   stored cell on the same flux surface;
3. continues from them (`hermes-3 -d <dir> restart ...`);
4. runs `summary.py`.

Measured (8 processes unless stated, machine shared):

| Case | Plasma time | Wall time | Result |
| --- | --- | --- | --- |
| Same grid, 8 processes, no change | 300 ms | 3 min 1 s | identical to the reference |
| Another GOAT grid of the same equilibrium (3,407 cells), 4 processes | 400 ms | 2 min 46 s | within 0.9 % of the reference |
| `power_core=20e6` | 400 ms | 6 min 6 s | steady: separatrix T_e 221 eV, outer-target peak 16.9 MW/m^2 |
| Grid of the public "VH" equilibrium (other strike points) | 32 ms | stopped at 20 min | not converged |

Process counts: any divisor of 48. `TIMEOUT` (seconds, default 3600) stops a run
that has gone wrong.

### 7c. Cold start (not under 10 minutes)

    mkdir -p $WORK/runs/cold && cp sparc-lmode/BOUT.inp $WORK/runs/cold/
    mpirun -np 8 $WORK/hermes-3-unstructured/build/hermes-3 -d $WORK/runs/cold \
        mesh:file=$PWD/sparc_lmode.grd.nc > $WORK/runs/cold/run.log 2>&1

starts from the uniform values in the input. It reached 5.5 ms of plasma time
after 10 minutes. The first milliseconds are the expensive part: 17 MW goes into
a cold plasma at rest, and until the flows and the divertor have settled the
solver takes steps of 5 to 30 microseconds. Starting from other uniform values,
from profiles in psi, or with charge exchange switched off at first does not
shorten it (tested). The stored state took 2.8 s of plasma time to reach,
through a detached phase. Use a cold start only for a device or set-up for which
no converged state exists.

### What a run prints

`BOUT.log.0` in the run directory has one line per time step,

    Time: 79600.6, timestep: 409.8, nl iter: 4, lin iter: 12, reason: 2

in units of 1/95788 ms (so `Time: 95788` is 1 ms), and one line per output with
the simulated time and the wall-clock time. When a step fails (Newton does not
converge in 7 iterations) the solver halves the step and tries again; the line
of the step that then succeeds ends with `SNES failures: <n>`. Thirty failures in
a row end the run. Up to one step in five needing a retry is normal in the first
milliseconds of a cold start.

The run is finished when the last line of `run.log` reports the run time and
`BOUT.dmp.*.nc` hold `nout + 1` time points. Do not open the output files of a
run that is still writing; wait for it or stop it.

### Continuing

    mpirun -np 8 $WORK/hermes-3-unstructured/build/hermes-3 -d <dir> restart \
        mesh:file=<the grid of that run> nout=20 timestep=957880

continues from the restart files in `<dir>` (10 ms between outputs here) and
overwrites the earlier output files; add `append` after `restart` to keep them.
Restart files are per process: continue with the same process count, or start a
new run from the old one with `STATE=<dir> STATE_GRID=<its grid> ./run_warm.sh`
(7b).

## 8. Results

    $WORK/venv/bin/python summary.py <run dir> sparc_lmode.grd.nc --json <run dir>/summary.json --figure <run dir>/summary.png

`summary.py` (read its header) prints and writes, in SI units with temperatures
in eV: how steady the run is, the power balance, the state on the core boundary
and on the outboard-midplane separatrix, the midplane profiles, and for each
divertor target the power, the ion flow and the peak heat flux with the plasma
state in front of it.

Steady state: `change_per_ms` gives the relative change per millisecond of the
electron content, the neutral content and the stored energy over the last output
interval. The stored state has all three below 2e-7; treat a run as steady when
all three are below 1e-5 (0.1 % per 100 ms). `history` in the JSON file has the
three totals at every output.

Reference values of the sample (the stored state on `sparc_lmode.grd.nc`):

| Quantity | Value |
| --- | --- |
| Plasma time of the stored state | 2.80 s |
| Outboard-midplane separatrix: n_e, T_e, T_i | 5.66e19 m^-3, 196 eV, 346 eV |
| Core boundary (psi_N 0.90): n_e, T_e, T_i, pressure | 8.0e19 m^-3, 2.42 keV, 2.77 keV, 66 kPa |
| Power: heating / sheath (electrons + ions) / neon / neutrals to walls | 17.0 / 4.96 + 5.26 / 2.68 / 0.45 MW |
| Remainder (hydrogen radiation, charge exchange) | 3.97 MW |
| Ion flow to the walls | 6.8e23 /s |
| Lower outer target: power; peak q on the surface; T_e there | 5.59 MW; 14.5 MW/m^2 (q_parallel 808 MW/m^2) at R 1.686, Z -1.339; 117 eV |
| Lower inner target: power; peak q on the surface; T_e there | 1.62 MW; 9.3 MW/m^2 at R 1.352, Z -1.166; 5.9 eV |
| Main-chamber wall (sheath) | 2.97 MW |

The outer target is attached and hot; the inner target is cold. The same numbers
are in `sparc-lmode/reference_summary.json`, which `--compare` reads.

For anything else, `tools/unstructured/usdata.py` reads every output variable
per cell in SI units:

    import sys; sys.path.insert(0, "<...>/tools/unstructured")
    import usdata
    run = usdata.Run("<run dir>", "sparc_lmode.grd.nc")
    Te = run("Te")            # [time, cell], eV
    R, Z = run.g["R"][:run.n], run.g["Z"][:run.n]   # cell centres, m

## 9. Another equilibrium

`make_grid.sh` takes any SPARC G-EQDSK file. What to check for a new file:

1. Units of psi. `geqdsk_to_rzpsi.py` tests them with Ampere's law (the current
   enclosed by the plasma boundary against the plasma current in the header)
   and prints which it found. It stops if the two disagree. A file without a
   plasma boundary needs `--psi-units` set by hand.
2. `RBTOR` printed by the script is written into GOAT's options automatically.
3. The vessel is GOAT's `Examples/SPARC/structure.dat`. The equilibrium box
   (`BOX` in `make_grid.sh`: R 1.1 to 2.5 m, Z -1.75 to 1.75 m) must contain it.
4. GOAT's log (`goat_case/goat.log`). "lines not allocated for cell" means psi
   lost precision near an X-point (do not reduce the digits written by
   `geqdsk_to_rzpsi.py`). Other GOAT failures are topology problems of the new
   equilibrium with the stored meshing options; they need a person.
5. The checks printed by `goat_to_bout.py`, as in section 5.

Then run as in 7b with the new grid: the stored state is mapped onto it along
flux surfaces. That is a good start only when the new equilibrium is close to the
sample's. For the public "VH" equilibrium (same plasma, other strike points) the
cells matched were 25 mm apart (median, against 0.1 mm between two grids of the
same equilibrium), and the run had reached only 32 ms after 20 minutes. Such a
case needs one long first run, wrapped in `timeout`; afterwards its own state
serves for fast restarts (`STATE=<that run> STATE_GRID=<its grid>`).

For a device other than SPARC you need its vessel as a GOAT `structure.dat`
and meshing options to match; start from the closest case in `goat/Examples`.

## 10. Problems seen, and what they mean

| Symptom | Cause and remedy |
| --- | --- |
| Hermes-3 stops at start with a list of unused options | BOUT++ refuses options nobody read: a misspelt name, or a top-level variable no expression uses. Fix the name or remove it |
| `Component type ... is not supported on flat-cell (unstructured) grids` | Only the components listed in `hermes-3.cxx` have been written for these grids. `sheath_boundary_unstructured` replaces the other sheath components |
| Load far above the process count, very slow run | `OPENBLAS_NUM_THREADS=1` was not set |
| Results differ after a rebuild | The build was configured without `-DHERMES_UPDATE_GIT_SUBMODULE=OFF` and dropped the BOUT++ patch, or with another slope limiter. `git -C external/BOUT-dev status --short` must list two modified files; the log prints the limiter |
| Error on restart about array sizes | The restart files belong to another process count or another grid. Start a new run from them with `STATE=<old dir> STATE_GRID=<old grid> ./run_warm.sh` (7a) |
| Steps stay at microseconds long after 3 ms, many `SNES failures` | Usually a plasma that is collapsing (too little power for the density and neon fraction). Look at `history` in `summary.json`; lower `neon_fraction` or `n_core` |
| `MatSolverType mumps` not available | PETSc was built without MUMPS (section 1) |
| HDF5 errors when reading output | The run was still writing. Read only finished runs |

## 11. What this model is not

- Neutrals are a fluid with a flux limiter, not kinetic. There are no molecules.
- Neon is a fixed fraction of the electron density with a coronal cooling curve.
  It is not transported and does not dilute the fuel.
- No drifts and no currents in this input (the branch has an experimental
  drift model; it is not used here).
- The pump and the gas valve are off by default, and their positions are
  placeholders. With them off, the particle content follows from `n_core`
  alone, and the separatrix density settles close to `n_core`. A fixed
  boundary density with `recycle < 1` or a pump lets particles, and the power
  they carry, cross the core boundary.
- Cross-field transport is the two numbers `D_perp` and `chi_perp`.
- The electron parallel heat flux has no kinetic limiter.
- The grid resolution was not varied for this case.
