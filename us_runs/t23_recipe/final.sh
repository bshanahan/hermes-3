#!/bin/bash
# Analysis and figures of the solver-recipe tests; copies the figures to paper/figures.
set -e
cd /home/shanahan/hermes-3-unstructured/us_runs
PY=/home/shanahan/hermes-3/sparc-dt-fast/venv/bin/python
D=t23_recipe/struct
G=t23_recipe/goat
$PY numerics/recipe_analyse.py numerics/recipe_struct.json $D/recipe $(ls -d $D/*)
$PY numerics/recipe_analyse.py numerics/recipe_goat.json $G/recipe $(ls -d $G/*)
$PY numerics/fig_recipe.py numerics/recipe_struct.json numerics/fig_recipe
$PY numerics/fig_grids.py numerics/recipe_struct.json numerics/recipe_goat.json numerics/fig_grids fastest
for f in fig_recipe fig_grids fig_convergence; do
  cp numerics/$f.pdf numerics/$f.png ../paper/figures/
done
echo "figures copied to paper/figures"
