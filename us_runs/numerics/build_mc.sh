#!/bin/bash
# Build a Hermes-3 binary with the MC slope limiter (Hermes-3's default), snapshot it, restore VanAlbada.
set -e
cd /home/shanahan/hermes-3-unstructured
echo "$(date +%T) configure MC"
nice -n 19 cmake -DHERMES_SLOPE_LIMITER=MC build > /dev/null
nice -n 19 cmake --build build --target hermes-3 -j 3
grep "Limiter=" build/include/hermes_build_config.hxx
mkdir -p us_runs/bin/mc
cp build/hermes-3 build/libhermes-3-lib.so us_runs/bin/mc/
cp us_runs/bin/v17/libbout++.so.5.2.0 us_runs/bin/mc/
echo "$(date +%T) snapshot in us_runs/bin/mc; restoring VanAlbada"
nice -n 19 cmake -DHERMES_SLOPE_LIMITER=VanAlbada build > /dev/null
nice -n 19 cmake --build build --target hermes-3 -j 3
grep "Limiter=" build/include/hermes_build_config.hxx
echo "$(date +%T) done"
