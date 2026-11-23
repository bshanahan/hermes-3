#!/bin/bash
# Where the GOAT grid's extra cost comes from, with the fastest recipe: skew correction on
# the plasma only, and none (two-point). Waits for PID $1.
cd /home/shanahan/hermes-3-unstructured/us_runs
if [ -n "$1" ]; then
  while kill -0 "$1" 2> /dev/null; do sleep 10; done
fi
F="solver:timestep_control=threshold_nonlinear_its"
bash t23_recipe/rqueue.sh \
  "t23_recipe/goat/fastest_plasma_only|t23_recipe/goat.inp|8|$F unstructured:skew_correction_neutrals=false|bin/v19" \
  "t23_recipe/goat/fastest_two_point|t23_recipe/goat.inp|8|$F unstructured:skew_correction=false|bin/v19"
