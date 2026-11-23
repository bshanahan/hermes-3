#!/bin/bash
# Baseline: the solver settings of Hermes-3's examples/tokamak-2D/recycling on the
# orthogonal TCV-X21 grid (8 processes), with van Albada (only the solver settings
# differ from the recipe) and with MC (Hermes-3 out of the box).
# Waits for the previous script (PID $1).
cd /home/shanahan/hermes-3-unstructured/us_runs
if [ -n "$1" ]; then
  while kill -0 "$1" 2> /dev/null; do sleep 10; done
fi
X="solver:max_nonlinear_iterations=10 solver:atol=1e-3 solver:maxf=2000 solver:matrix_free_operator=false"
X="$X solver:lag_jacobian=3 solver:jacobian_reuse=0 solver:max_snes_failures=20 solver:target_its=100"
X="$X solver:kP=0.233 solver:kI=0.133 petsc:snes_ksp_ew=false"
bash t23_recipe/rqueue.sh \
  "t23_recipe/struct/example_settings|t23_recipe/struct8.inp|8|$X|bin/v19" \
  "t23_recipe/struct/example_settings_mc|t23_recipe/struct8.inp|8|$X|bin/mc"
