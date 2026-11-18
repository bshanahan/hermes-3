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
