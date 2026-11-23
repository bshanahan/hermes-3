# TCV-X21 on a GOAT grid

Hermes-3 with flat-cell operators on a GOAT grid of the TCV-X21 reference
equilibrium (TCV shot 65402), compared with the public TCV-X21 dataset and
with the SOLPS-ITER case of Wang et al., Nucl. Fusion 64 (2024) 056040.
No drifts, deuterium only, fluid atoms, no carbon.

## Steps

1. Download the dataset (CC-BY 4.0) from https://github.com/SPCData/TCV-X21
   into `tcvx21/data/`. The files used are the eqdsk of the reference
   scenario and `TCV_forward_field.nc` / `TCV_reversed_field.nc`.
2. Make the GOAT input: `python goat/make_rzpsi.py <eqdsk> rzpsi.dat 534 1521`.
   The eqdsk stores psi in Wb, not Wb/rad; the script detects this with
   Ampere's law.
3. Run GOAT from a copy of its `Examples/TCV` case with three options changed:
   `goat.magneticfieldfilepath 'rzpsi.dat'`, `goat.mf.readmeth 'rzpsi'`,
   `goat.mf.RBtor '-0.8204806'`.
4. Convert: `python ../tools/unstructured/goat_to_bout.py traduit.out.b2us grids/tcvx21_goat.grd.nc`.
5. Run from `us_runs/`: `mpirun -np 8 ../build/hermes-3 -d t12_tcvx21`
   with `BOUT.inp.stage1` copied to `BOUT.inp` (about 90 s).
6. Compare: `python tcvx21/compare.py fig.png "label=t12_tcvx21"` and
   `python tcvx21/analyse.py t12_tcvx21`.

`structured/make_grid.py` makes an orthogonal hypnotoad grid of the same
equilibrium and wall for comparison.

## Inputs in `t12_tcvx21/`

| Input | Fuelling | Wall |
| --- | --- | --- |
| `BOUT.inp.stage1` | density 1.2e19 m^-3 held at the core boundary | absorbs 1% of ions and neutrals |
| `BOUT.inp.stage2` | puff 1.36e20 atoms/s (the SOLPS-ITER rate), no core flux | as above; the plasma is not sustained |
| `puff_9e20/BOUT.inp` | puff 9e20 atoms/s | as above |
| `ions_only_1.36e20/BOUT.inp` | puff 1.36e20 atoms/s | absorbs 1% of ions, reflects neutrals |

The valve position and the reading of the paper's recycling coefficient
(R = 0.99) are assumptions; the paper does not give them in full.

## Orthogonal grid

`t13_tcvx21_struct/BOUT.inp` is the stage-1 case on the hypnotoad grid, run
with the same binary (`mpirun -np 13`, so that the branch cuts fall on
processor boundaries). Compare both grids with

    python tcvx21/compare.py fig.png "GOAT grid=t12_tcvx21" \
        "orthogonal grid=t13_tcvx21_struct::tcvx21/structured/tcvx21_32x104.grd.nc"

and get the balances with `python tcvx21/balance_struct.py t13_tcvx21_struct`.
The two cases are not identical at the boundaries: the structured grid ends on
flux surfaces (psiN 1.24 and 0.75) instead of the wall, its target rows are
perpendicular to the flux surfaces, and its targets do not absorb neutrals.

## Further variants and checks

| Input | What it changes |
| --- | --- |
| `t12_tcvx21/wall_decay_3cm/BOUT.inp` | stage 1 with `wall_decay_length = 0.03` for ions and electrons |
| `t15_tcvx21_carbon/BOUT.inp` | stage 1 with carbon sputtered from every wall face (`wall_sputtering`, yield 3.5% of the deuterium ion flow), charge states up to C4+ |

The carbon case reads OpenADAS rates from `json_database/` in the directory
Hermes-3 is started from: link it there (`ln -s ../json_database us_runs/`).
It restarts from a stage-1 state to which trace carbon has been added.
`python tcvx21/carbon.py <run dir>` prints the carbon content, source and
radiation (rebuilt from the same tables; the ADAS reactions write no
radiation diagnostic).

Sputtering by the atoms that reach the wall as well (`from = d+, d`) gives
several times more carbon than SOLPS-ITER has, because fluid atoms touch the
wall about ten times per ionisation.

`python tcvx21/reference.py check <run dir> tcvx21/reference/stage1_default.json`
compares a stage-1 run with the recorded numbers: separatrix values, flows and
the distance of each profile to the measurements. `stage1_two_point.json` is
the record for `unstructured:skew_correction=false`.

## Carbon variants

| Input | What it changes from `t15_tcvx21_carbon` |
| --- | --- |
| `t16_tcvx21_carbon_Ti/BOUT.inp` | every carbon ion evolves its own pressure (`evolve_pressure`) instead of taking the deuterium ion temperature |
| `t17_tcvx21_carbon_atom_momentum/BOUT.inp` | carbon atoms evolve their parallel momentum; `recycling:pump_momentum_sink = true` so that atoms absorbed by the wall take their momentum with them |

Both restart from the state of `t15_tcvx21_carbon`. For `t16` the restart files
need the carbon ion pressures (`P = n T_D+`).
`python tcvx21/carbon_variants.py <label>=<run dir> ...` prints them side by side.

## Drifts

`t18_tcvx21_drifts/BOUT.inp` is stage 1 with the component `drifts_unstructured`
(see `include/drifts_unstructured.hxx`): E x B and magnetic drifts of every
charged species through the cell faces, and a potential found from current
continuity (Ohm's law along the field, sheath currents at the wall, the
current of the magnetic drifts, and a small anomalous perpendicular current).
It needs a grid file with the drift geometry (`us_face_drift`, `us_face_curl`
and vertex stencils on wall faces): `grids/tcvx21_goat_drifts.grd.nc` is the
standard grid written again with the current `goat_to_bout.py`.

Steps:

1. `python tcvx21/add_phi.py <dir>` adds a first guess of the potential
   (3 T_e) to the restart files of a stage-1 run.
2. Run with `drifts:exb_advection=false drifts:magnetic_drift=false` until
   steady: only the currents act (`potential_only`).
3. Restart from that state with the drifts on. `unstructured:field_sign = 1` is
   the field of the grid file, in which the ion grad-B drift points to the
   X-point (the dataset's forward field); `-1` is the reversed field.

`python tcvx21/drifts.py <run dir> [<reference run>] [figure.png]` prints the
potential, the current balance and the drift flows.

Checks that can fail:

* `tests/unstructured/test_operators.py`: the drift geometry of the grid file
  against a quadrature of the analytic field (sign and size), convergence of
  the E x B flow, conservation.
* With `drifts:conductivity_factor=1e-12 drifts:perp_mobility=1e-14` and both
  drifts off, every sheath floats and the run must stay on the stage-1
  solution (it does, to 0.01%). With the real conductivity it does not: a
  thermoelectric current of about 220 A flows between the targets.

## Resolution scan

`t19_tcvx21_resolution/<grid>/{nodrift,forward}` repeat stage 1 and the
forward-field drift case on GOAT grids made with other cell sizes
(`~/goat/Runs/TCV_X21_<grid>`): `coarse` and `fine` have half and twice the
cells along each flux tube. `radial` has twice the flux tubes across the
midplane (GOAT option `gg.radref.LB.Lmaxinf = 0.008`), but its thinnest cells
are twenty times smaller than those of the standard grid and neither case has
run on it yet: the implicit solver stays at steps of about 0.1 microseconds,
failing on the neutral equations.
A state is moved between grids with `tools/unstructured/regrid_restart.py`.
`python tcvx21/scan.py <label>=<run dir>::<grid file> ...` prints the runs
side by side; it picks one cell per flux tube along each measurement line,
which works at any resolution.
