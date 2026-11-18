# Flat-cell ("unstructured") grids for Hermes-3

Experimental. Cells of any shape are stored along the BOUT++ y index
(`nx = 1`, `MXG = 0`, `MZ = 1`), with one ghost cell per boundary face that
holds the value at that face. The operators are finite-volume and work on the
face list in the grid file; see `include/unstructured.hxx` and the header of
`flatgrid.py` for the grid-file contents.

| File | Purpose |
| --- | --- |
| `flatgrid.py` | writes a flat-cell grid file from cells, faces and a psi map |
| `goat_to_bout.py` | converts a GOAT grid (`traduit.out.b2us`) to a flat-cell grid |
| `b2us.py` | reader for `traduit.out.b2us` |
| `hypnotoad_to_flat.py` | recasts a structured hypnotoad grid as a flat-cell grid with the same cells, for verification |
| `compare_same_cells.py` | compares a structured run with a flat-cell run on the same cells |
| `usdata.py` | reads flat-cell run output in SI units |
| `usops.py` | numpy reference implementation of the operators, used to check the C++ |
| `BOUT-dev.patch` | two changes to BOUT-dev that the flat-cell runs need |

## BOUT-dev patch

The submodule pointer is unchanged. Apply the patch before building:

    git -C external/BOUT-dev apply ../../tools/unstructured/BOUT-dev.patch

It touches `src/solver/impls/snes/snes.cxx` and
`src/solver/petsc_preconditioner.cxx`.

## Input

A flat-cell run uses the ordinary components plus
`sheath_boundary_unstructured` and, for gas valves, `wall_puff`.
`us_runs/t12_tcvx21/BOUT.inp.stage1` is a complete example.
