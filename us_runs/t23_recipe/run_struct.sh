#!/bin/bash
# Structured-grid recipe tests: each job changes ONE thing from the recipe (struct8.inp).
cd /home/shanahan/hermes-3-unstructured/us_runs
I=t23_recipe/struct8.inp
D=t23_recipe/struct
REF=LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/blas
bash t23_recipe/rqueue.sh \
  "$D/recipe|$I|8||" \
  "$D/jac_every_step|$I|8|solver:jacobian_reuse=0|" \
  "$D/lim_MC|$I|8||bin/mc" \
  "$D/lim_MinMod|$I|8||bin/lim_MinMod" \
  "$D/lim_Superbee|$I|8||bin/lim_Superbee" \
  "$D/lim_Upwind|$I|8||bin/lim_Upwind" \
  "$D/lim_WENO3|$I|8||bin/lim_WENO3" \
  "$D/coulomb_orig|$I|8||bin/coulomb_orig" \
  "$D/no_ew|$I|8|petsc:snes_ksp_ew=false|" \
  "$D/predictor_const|$I|8|solver:predictor:type=constant|" \
  "$D/threshold_dt|$I|8|solver:timestep_control=threshold_nonlinear_its|" \
  "$D/atol_1e-3|$I|8|solver:atol=1e-3|" \
  "$D/ref_blas|$I|8|||$REF" \
  "$D/jac_persists|$I|8|solver:jacobian_reuse=0 solver:jacobian_persists=true|" \
  "$D/jac_every_newton|$I|8|solver:jacobian_reuse=0 solver:lag_jacobian=1|" \
  "$D/no_limiter_lag|$I|8|d:lag_limiter_gradient=false|" \
  "$D/fd_operator|$I|8|solver:matrix_free_operator=false|" \
  "$D/bjacobi_ilu|$I|8|solver:pc_type=bjacobi|" \
  "$D/recipe_rep|$I|8||"
