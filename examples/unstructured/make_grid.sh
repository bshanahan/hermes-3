#!/bin/bash
# Flat-cell (unstructured) Hermes-3 grid of SPARC from a G-EQDSK equilibrium, with GOAT.
#
# usage: make_grid.sh <GOAT directory> <equilibrium.geqdsk> <output.grd.nc> [work directory]
#
#   <GOAT directory>   a built checkout of https://github.com/TME-TFE/goat
#                      (executables/goat.exe must exist)
#   [work directory]   where GOAT runs; default ./goat_case (created, files overwritten)
#
# Environment:
#   PYTHON        python with numpy, scipy and netCDF4 (default: python3)
#   GOAT_THREADS  OpenMP threads for GOAT (default: 4)
#
# Steps:
#   1. geqdsk_to_rzpsi.py writes GOAT's equilibrium file rzpsi.dat (psi in Wb,
#      separatrix at zero) on the box of GOAT's SPARC example, and reports R*Btor.
#   2. The vessel (structure.dat) and the meshing options (GOAToptions.dat) are those
#      of GOAT's own SPARC example; only 'goat.mf.RBtor' is set from the equilibrium.
#   3. goat.exe writes the mesh, traduit.out.b2us.
#   4. goat_to_bout.py converts it to the grid file read by Hermes-3 and checks it.
set -euo pipefail

if [ $# -lt 3 ]; then
  sed -n '2,20p' "$0"; exit 2
fi
GOAT=$(cd "$1" && pwd)
GEQDSK=$(realpath "$2")
OUT=$(realpath -m "$3")
WORK=$(realpath -m "${4:-goat_case}")
PYTHON=${PYTHON:-python3}
TOOLS=$(cd "$(dirname "$0")/../../tools/unstructured" && pwd)

# R, Z box of the equilibrium given to GOAT [m]: contains the SPARC vessel, leaves out the coils
BOX="1.1 2.5 -1.75 1.75"
POINTS="1025 1025"

[ -x "$GOAT/executables/goat.exe" ] || { echo "error: $GOAT/executables/goat.exe not found: build GOAT first" >&2; exit 1; }
[ -f "$GOAT/Examples/SPARC/structure.dat" ] || { echo "error: $GOAT/Examples/SPARC not found" >&2; exit 1; }
command -v tcsh > /dev/null || { echo "error: tcsh is needed to run GOAT" >&2; exit 1; }

mkdir -p "$WORK"
cp "$GOAT/Examples/SPARC/structure.dat" "$GOAT/Examples/SPARC/GOAToptions.dat" "$WORK/"

echo "== 1. equilibrium -> $WORK/rzpsi.dat"
"$PYTHON" "$TOOLS/geqdsk_to_rzpsi.py" "$GEQDSK" "$WORK/rzpsi.dat" --box $BOX --points $POINTS | tee "$WORK/geqdsk_to_rzpsi.log"
RBTOR=$(grep '^RBTOR=' "$WORK/geqdsk_to_rzpsi.log" | cut -d= -f2)
[ -n "$RBTOR" ] || { echo "error: no RBTOR line from geqdsk_to_rzpsi.py" >&2; exit 1; }

echo "== 2. GOAT options: goat.mf.RBtor = $RBTOR"
grep -q "^'goat.mf.RBtor'" "$WORK/GOAToptions.dat" || { echo "error: no goat.mf.RBtor line in GOAToptions.dat" >&2; exit 1; }
sed -i "s|^'goat.mf.RBtor'.*|'goat.mf.RBtor'                '$RBTOR'|" "$WORK/GOAToptions.dat"

echo "== 3. GOAT (log: $WORK/goat.log)"
rm -f "$WORK/traduit.out.b2us"
tcsh -c "cd '$GOAT'; source setup.csh > /dev/null; setenv OMP_NUM_THREADS ${GOAT_THREADS:-4}; cd '$WORK'; goat.exe" > "$WORK/goat.log" 2>&1 || true
[ -s "$WORK/traduit.out.b2us" ] || { echo "error: GOAT wrote no traduit.out.b2us; see $WORK/goat.log" >&2; tail -20 "$WORK/goat.log" >&2; exit 1; }

echo "== 4. GOAT mesh -> $OUT"
"$PYTHON" "$TOOLS/goat_to_bout.py" "$WORK/traduit.out.b2us" "$OUT"
