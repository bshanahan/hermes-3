#!/bin/bash
# Binaries for the solver-recipe tests (t23_recipe): one per slope limiter, and one
# with the original (discontinuous) e-i Coulomb logarithm. Each is snapshotted in
# us_runs/bin/<name>; afterwards build/ is restored to the production state (= bin/v17).
set -e
cd /home/shanahan/hermes-3-unstructured
J=12
snap() {
  mkdir -p us_runs/bin/$1
  cp build/hermes-3 build/libhermes-3-lib.so us_runs/bin/$1/
  cp us_runs/bin/v17/libbout++.so.5.2.0 us_runs/bin/$1/
  md5sum us_runs/bin/$1/hermes-3 us_runs/bin/$1/libhermes-3-lib.so
}
for L in MinMod Superbee Upwind WENO3; do
  echo "$(date +%T) limiter $L"
  cmake -DHERMES_SLOPE_LIMITER=$L build > /dev/null
  cmake --build build --target hermes-3 -j $J > us_runs/t23_recipe/build_$L.log 2>&1
  grep "Limiter=" build/include/hermes_build_config.hxx
  snap lim_$L
done
echo "$(date +%T) VanAlbada + original Coulomb logarithm"
cmake -DHERMES_SLOPE_LIMITER=VanAlbada build > /dev/null
cp src/braginskii_collisions.cxx us_runs/t23_recipe/braginskii_collisions.cxx.current
git show 5aab5b8a^:src/braginskii_collisions.cxx > src/braginskii_collisions.cxx
cmake --build build --target hermes-3 -j $J > us_runs/t23_recipe/build_coulomb.log 2>&1
grep "Limiter=" build/include/hermes_build_config.hxx
snap coulomb_orig
cp us_runs/t23_recipe/braginskii_collisions.cxx.current src/braginskii_collisions.cxx
echo "$(date +%T) restore the production build"
cmake --build build --target hermes-3 -j $J > us_runs/t23_recipe/build_restore.log 2>&1
grep "Limiter=" build/include/hermes_build_config.hxx
git diff --stat -- src/braginskii_collisions.cxx
md5sum build/hermes-3 build/libhermes-3-lib.so us_runs/bin/v17/hermes-3 us_runs/bin/v17/libhermes-3-lib.so
echo "$(date +%T) done"
