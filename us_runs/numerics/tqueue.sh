#!/bin/bash
# Timing queue (run from us_runs/). job = dir|input|nranks|extra|restart_from
# Runs one job at a time; logs load before and after each run to tqueue.log in the first job's parent.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
BIN=${HERMES_BIN:-bin/v17}
for job in "$@"; do
  IFS='|' read -r dir input np extra from <<< "$job"
  mkdir -p "$dir"
  cp "$input" "$dir/BOUT.inp"
  rm -f "$dir"/BOUT.dmp.*.nc "$dir"/BOUT.log.* 
  rs=""
  if [ -n "$from" ]; then cp "$from"/BOUT.restart.*.nc "$dir/"; rs="restart"; fi
  log="$(dirname "$dir")/tqueue.log"
  echo "$(date '+%F %T') start $dir np=$np extra=[$extra] restart=[$from] load: $(cut -d' ' -f1-3 /proc/loadavg)" >> "$log"
  t0=$(date +%s.%N)
  timeout 2h mpirun -np "$np" "$BIN/hermes-3" -d "$dir" $rs $extra > "$dir/run.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  echo "$(date '+%F %T') end   $dir exit=$rc wall=$(echo "$t1 - $t0" | bc) s failures=$(grep -a -c 'reason: -' "$dir/BOUT.log.0") steps=$(grep -a -c 'reason: [0-9]' "$dir/BOUT.log.0") $(grep -a 'Run time' "$dir/run.log" | tail -1) load: $(cut -d' ' -f1-3 /proc/loadavg)" >> "$log"
done
