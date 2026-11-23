#!/bin/bash
# GOAT-grid runs with the solver recipe (binary v19 = new bound on the skew correction,
# with the field's scale by kind: harmonic mean / 1 for logarithms / rms for signed fields).
# Waits for the structured queue (PID given as $1) so that no two timed runs overlap.
cd /home/shanahan/hermes-3-unstructured/us_runs
if [ -n "$1" ]; then
  while kill -0 "$1" 2> /dev/null; do sleep 10; done
fi
G=t23_recipe/goat.inp
D=t23_recipe/goat
OLD="unstructured:skew_smooth=0.05 unstructured:skew_harmonic_mean=false"
bash t23_recipe/rqueue.sh \
  "$D/recipe|$G|8||bin/v19" \
  "$D/old_bound|$G|8|$OLD|bin/v19" \
  "$D/plasma_only|$G|8|unstructured:skew_correction_neutrals=false|bin/v19" \
  "$D/two_point|$G|8|unstructured:skew_correction=false|bin/v19" \
  "$D/jac_every_step|$G|8|solver:jacobian_reuse=0|bin/v19" \
  "$D/recipe_rep|$G|8||bin/v19" \
  "t23_recipe/struct/recipe_v19|t23_recipe/struct8.inp|8||bin/v19"
