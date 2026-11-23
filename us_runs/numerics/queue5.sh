#!/bin/bash
# after queue 4: MC-limiter cold start, then PVODE (variable-order BDF) for at most 10 minutes each
cd /home/shanahan/hermes-3-unstructured/us_runs
while kill -0 2729710 2>/dev/null; do sleep 15; done
until [ -x bin/mc/hermes-3 ] && grep -q "done" t22_timing/build_mc.log; do sleep 30; done
HERMES_BIN=bin/mc /tmp/claude-1000/-home-shanahan-spielplatz-wsl--venv/a5aac4d2-4fa1-4f58-9c4c-435a93df0aa2/scratchpad/paper_tests/tqueue.sh 't22_timing/flat_default_mc|t22_timing/flat.inp|8||'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for case in cold restart; do
  d=t22_timing/pvode_$case; mkdir -p $d; rm -f $d/BOUT.dmp.* $d/BOUT.log.*; cp t22_timing/flat.inp $d/BOUT.inp
  rs=""; [ $case = restart ] && { cp t22_timing/restart_np8/BOUT.restart.*.nc $d/; rs=restart; }
  echo "$(date '+%F %T') start $d" >> t22_timing/tqueue.log
  t0=$(date +%s)
  timeout 600 mpirun -np 8 bin/v17/hermes-3 -d $d $rs solver:type=pvode solver:mxstep=100000000 input:error_on_unused_options=false > $d/run.log 2>&1
  echo "$(date '+%F %T') end   $d exit=$? wall=$(( $(date +%s) - t0 )) s" >> t22_timing/tqueue.log
done
