#!/bin/bash
# The fastest recipe found on the orthogonal grid (the recipe with threshold time-step
# control instead of PID) on both grids, twice each, binary v19.
# Waits for the previous script (PID $1) so that no two timed runs overlap.
cd /home/shanahan/hermes-3-unstructured/us_runs
if [ -n "$1" ]; then
  while kill -0 "$1" 2> /dev/null; do sleep 10; done
fi
F="solver:timestep_control=threshold_nonlinear_its"
bash t23_recipe/rqueue.sh \
  "t23_recipe/struct/fastest|t23_recipe/struct8.inp|8|$F|bin/v19" \
  "t23_recipe/goat/fastest|t23_recipe/goat.inp|8|$F|bin/v19" \
  "t23_recipe/struct/fastest_rep|t23_recipe/struct8.inp|8|$F|bin/v19" \
  "t23_recipe/goat/fastest_rep|t23_recipe/goat.inp|8|$F|bin/v19"
