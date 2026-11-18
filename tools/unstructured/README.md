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
| `annulus.py` | analytic test grid: an annulus of circular flux surfaces, optionally skewed |
| `skew_report.py` | size of the cross-field flux error on skewed cells of a GOAT grid, with and without the correction |
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

## Cross-field operators on skewed cells

The two-point flow through a face uses the two cell centres only. It is exact
when the line joining them is normal to the face. Otherwise the gradient along
the face is missed: on the TCV-X21 GOAT grid the error for a gradient along the
face exceeds 25% of `D |grad f| A` on a tenth of the interior faces.

Grid files written by `flatgrid.py` therefore carry, for every face, the
coefficient of a vertex term and the cells from which the two vertex values
are interpolated (a least-squares linear fit). With

    [unstructured]
    skew_correction = true            # default if the grid file has the data
    skew_correction_neutrals = true   # default: same as skew_correction
    skew_limit = 1.0                  # correction <= this multiple of the two-point term
    skew_smooth = 0.05                # smoothing of that bound

the cross-field operators of the plasma (`anomalous_diffusion`) and of the
neutrals (`neutral_mixed`) add that term. A linear field then has exact flows,
and the flow error on a skewed grid falls at second order
(`tests/unstructured`).

The correction is bounded by `skew_limit` times the two-point term, so that
the flow between two cells cannot run strongly from the lower to the higher
value; without a bound the neutral density goes negative on the TCV-X21 case.
The bound is smoothed: `skew_smooth` times the mean size of the field counts
as a difference between the two cells. With `skew_smooth = 0` the bound has a
kink wherever two cells hold the same value and the implicit solver takes very
small steps once the neutrals are corrected.

## Cross-field loss to the wall

By default no cross-field flux crosses a wall face. `wall_decay_length` [m] in
a species section (next to `anomalous_D`) adds the loss `D n / length` for
particles and `chi n T / length` for heat through the part of each wall face
that faces across the flux surfaces. The ions lost are recycled by the
`recycling` component like those that arrive through the sheath.

## Supported components

On a flat-cell grid Hermes-3 stops if the input names a component that has not
been written or checked for these grids (the list is in `hermes-3.cxx`).
`unstructured:allow_untested_components = true` turns the error into a warning.

## Tests

    python tests/unstructured/test_operators.py        # grid tools and reference operators
    tests/integrated/unstructured-diffusion/runtest    # Hermes-3 against the reference, 1 and 3 processes

The integrated test is registered with CTest. Set `HERMES_EXE` to run it by
hand from the source tree.
