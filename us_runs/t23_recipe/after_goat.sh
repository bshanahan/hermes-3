#!/bin/bash
# After the GOAT queue (PID $1): agreement test with v19, and the golden record for the new bound.
cd /home/shanahan/hermes-3-unstructured/us_runs
while kill -0 "$1" 2> /dev/null; do sleep 10; done
PY="/home/shanahan/hermes-3/sparc-dt-fast/venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "== $(date +%T) agreement test, v19"
$PY numerics/agreement.py setup v19 && $PY numerics/agreement.py analyse v19
REF=tcvx21/reference
echo "== $(date +%T) v19 with the first bound against the record made with v17"
$PY -I tcvx21/reference.py check t23_recipe/goat/old_bound $REF/fixed_160kW_nodrift.json 1e-3
echo "== $(date +%T) new bound against the record made with v17 (expected to differ)"
$PY -I tcvx21/reference.py check t23_recipe/goat/recipe $REF/fixed_160kW_nodrift.json 1e-3
echo "== $(date +%T) record for the new bound"
$PY -I tcvx21/reference.py record t23_recipe/goat/recipe $REF/fixed_160kW_nodrift_v19.json \
  "v19 (bound on the skew correction: smooth 1, per-field scale), solver recipe, 8 ranks, cold start to 300 ms"
echo "== $(date +%T) repeat run against the new record"
$PY -I tcvx21/reference.py check t23_recipe/goat/recipe_rep $REF/fixed_160kW_nodrift_v19.json 1e-3
echo "== $(date +%T) two-point against the new record (effect of the skew correction)"
$PY -I tcvx21/reference.py check t23_recipe/goat/two_point $REF/fixed_160kW_nodrift_v19.json 1e-3
echo "== $(date +%T) done"
