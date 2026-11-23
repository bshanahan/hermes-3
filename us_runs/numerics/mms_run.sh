#!/bin/bash
# usage (from us_runs/numerics/mms): ../mms_run.sh <nranks> <tag> [<tag> ...]
# Runs every case <tag> (an .inp written by mms_setup.py) with the three operator settings.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
np=$1; shift
for tag in "$@"; do
  for cfg in default unbounded two_point; do
    case $cfg in
      default)   extra="" ;;
      unbounded) extra="unstructured:skew_limit=1e12" ;;
      two_point) extra="unstructured:skew_correction=false" ;;
    esac
    dir="${tag}_${cfg}"
    if grep -q "Run time" "$dir/run.log" 2>/dev/null; then continue; fi
    mkdir -p "$dir"; cp "$tag.inp" "$dir/BOUT.inp"; rm -f "$dir"/BOUT.dmp.*
    t0=$(date +%s.%N)
    timeout 3600 mpirun -np "$np" ../../bin/v17/hermes-3 -d "$dir" $extra > "$dir/run.log" 2>&1
    rc=$?
    echo "$(date +%T) $dir np=$np exit=$rc wall=$(echo "$(date +%s.%N) - $t0" | bc | cut -c1-6)s failures=$(grep -a -c 'reason: -' $dir/BOUT.log.0)"
  done
done
