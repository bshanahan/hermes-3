#!/bin/bash
# Clean rerun of the block-Jacobi ILU(0) test (the first one shared the machine with
# diagnostic runs). 40 ms is enough: the steady state is reached at ~8 ms; the time to
# steady state is measured against the recipe's final state. Waits for PID $1.
cd /home/shanahan/hermes-3-unstructured/us_runs
if [ -n "$1" ]; then
  while kill -0 "$1" 2> /dev/null; do sleep 10; done
fi
bash t23_recipe/rqueue.sh \
  "t23_recipe/struct/bjacobi_ilu_clean|t23_recipe/struct8.inp|8|solver:pc_type=bjacobi nout=40|bin/v19"
