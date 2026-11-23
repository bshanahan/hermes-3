#!/bin/bash
# usage (from us_runs/numerics/mms): ../mms_run2.sh <nranks> <binary dir> "<cfg> <cfg> ..." <tag> [<tag> ...]
# Runs every case <tag> (an .inp written by mms_setup.py) with the named bound settings;
# output in <tag>_<cfg>. Settings (binary v18 or later):
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
np=$1; bin=$2; cfgs=$3; shift 3
for tag in "$@"; do
  for cfg in $cfgs; do
    case $cfg in
      new)       extra="" ;;                                       # code default: smooth 1, harmonic, p 2
      k05)       extra="unstructured:skew_smooth=0.5" ;;
      k025)      extra="unstructured:skew_smooth=0.25" ;;
      k1p8)      extra="unstructured:skew_bound_power=8" ;;
      old)       extra="unstructured:skew_smooth=0.05 unstructured:skew_harmonic_mean=false" ;;
      unbounded) extra="unstructured:skew_limit=1e12" ;;
      two_point) extra="unstructured:skew_correction=false" ;;
      *) echo "unknown cfg $cfg"; exit 1 ;;
    esac
    dir="${tag}_${cfg}"
    if grep -q "Run time" "$dir/run.log" 2>/dev/null; then continue; fi
    mkdir -p "$dir"; cp "$tag.inp" "$dir/BOUT.inp"; rm -f "$dir"/BOUT.dmp.*
    t0=$(date +%s.%N)
    timeout 3600 mpirun -np "$np" ../../$bin/hermes-3 -d "$dir" $extra > "$dir/run.log" 2>&1
    rc=$?
    echo "$(date +%T) $dir np=$np exit=$rc wall=$(echo "$(date +%s.%N) - $t0" | bc | cut -c1-6)s failures=$(grep -a -c 'reason: -' $dir/BOUT.log.0)"
  done
done
