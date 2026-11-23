#!/bin/bash
# Solver-recipe test queue (run from us_runs/).
#   job = dir|input|nranks|extra options|binary dir (default bin/v17)|environment (VAR=value ...)
# One job at a time, each a cold start, stopped after 1 hour (exit 124).
# Start and end lines, with the machine load, go to t23_recipe/rqueue.log.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
LOG=t23_recipe/rqueue.log
for job in "$@"; do
  IFS='|' read -r dir input np extra bin envs <<< "$job"
  bin=${bin:-bin/v17}
  mkdir -p "$dir"
  cp "$input" "$dir/BOUT.inp"
  rm -f "$dir"/BOUT.dmp.*.nc "$dir"/BOUT.log.* "$dir"/BOUT.restart.*.nc
  echo "$(date '+%F %T') start $dir np=$np bin=$bin extra=[$extra] env=[$envs] load: $(cut -d' ' -f1-3 /proc/loadavg)" >> $LOG
  t0=$(date +%s.%N)
  env $envs timeout 3600 mpirun -np "$np" "$bin/hermes-3" -d "$dir" $extra > "$dir/run.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  echo "$(date '+%F %T') end   $dir exit=$rc wall=$(echo "$t1 - $t0" | bc) s failures=$(grep -a -c 'reason: -' "$dir/BOUT.log.0") steps=$(grep -a -c 'reason: [0-9]' "$dir/BOUT.log.0") $(grep -a 'Run time' "$dir/run.log" | tail -1) load: $(cut -d' ' -f1-3 /proc/loadavg)" >> $LOG
done
echo "$(date '+%F %T') queue finished" >> $LOG
