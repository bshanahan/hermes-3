#!/bin/bash
# Run the SPARC sample from its converged state (restart files in sparc-lmode/) on any
# grid and process count.
#
# usage: run_warm.sh <run directory> [Hermes-3 options ...]
#
# Options after the run directory go to Hermes-3, for example
#     run_warm.sh run1 power_core=12e6 n_core=8e19
#     run_warm.sh run2 nout=10 timestep=957880          (10 outputs, 10 ms apart)
#
# Environment (defaults in brackets):
#   HERMES   Hermes-3 executable            [<this repository>/build/hermes-3]
#   PYTHON   python with numpy, scipy, netCDF4 [python3]
#   NP       MPI processes, a divisor of 48  [8]
#   GRID     grid file                      [sparc_lmode.grd.nc next to this script]
#   INPUT    input file                     [sparc-lmode/BOUT.inp next to this script]
#   STATE    directory with the restart files of the stored state [sparc-lmode next to this script]
#   STATE_GRID  grid of the stored state     [sparc_lmode.grd.nc next to this script]
#   MPIRUN   MPI launcher                   [mpirun]
#   TIMEOUT  stop the run after this many seconds [3600]
#
# Steps:
#   1. copy the input into the run directory and run one tiny step, only to create
#      restart files for this grid and process count (init.log);
#   2. fill them with the stored state, mapped onto this grid along flux surfaces
#      (tools/unstructured/regrid_restart.py; regrid.log);
#   3. continue from them with the options given (run.log);
#   4. write <run directory>/summary.json with summary.py and print the main numbers.
# The run directory must not exist yet, or be empty.
set -euo pipefail

if [ $# -lt 1 ]; then
  sed -n '2,28p' "$0"; exit 2
fi
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
RUN=$(realpath -m "$1"); shift
HERMES=${HERMES:-$REPO/build/hermes-3}
PYTHON=${PYTHON:-python3}
NP=${NP:-8}
GRID=$(realpath "${GRID:-$HERE/sparc_lmode.grd.nc}")
INPUT=$(realpath "${INPUT:-$HERE/sparc-lmode/BOUT.inp}")
STATE=$(realpath "${STATE:-$HERE/sparc-lmode}")
STATE_GRID=$(realpath "${STATE_GRID:-$HERE/sparc_lmode.grd.nc}")
MPIRUN=${MPIRUN:-mpirun}
TIMEOUT=${TIMEOUT:-3600}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

[ -x "$HERMES" ] || { echo "error: Hermes-3 executable $HERMES not found (set HERMES)" >&2; exit 1; }
ls "$STATE"/BOUT.restart.*.nc > /dev/null 2>&1 || { echo "error: no restart files in $STATE" >&2; exit 1; }
[ -f "$STATE_GRID" ] || { echo "error: $STATE_GRID not found" >&2; exit 1; }
if [ -e "$RUN" ] && [ -n "$(ls -A "$RUN")" ]; then
  echo "error: $RUN exists and is not empty" >&2; exit 1
fi
mkdir -p "$RUN"
cp "$INPUT" "$RUN/BOUT.inp"

# Options for step 1: everything given, except the output settings
INIT_OPTS=()
for o in "$@"; do
  case "$o" in nout=*|timestep=*) ;; *) INIT_OPTS+=("$o") ;; esac
done

echo "== 1. restart files for $(basename "$GRID") on $NP processes"
"$MPIRUN" -np "$NP" "$HERMES" -d "$RUN" mesh:file="$GRID" ${INIT_OPTS[@]+"${INIT_OPTS[@]}"} nout=1 timestep=1e-3 \
  > "$RUN/init.log" 2>&1 || { echo "error: see $RUN/init.log" >&2; tail -20 "$RUN/init.log" >&2; exit 1; }
rm -f "$RUN"/BOUT.dmp.*.nc

echo "== 2. stored state -> restart files"
"$PYTHON" "$REPO/tools/unstructured/regrid_restart.py" "$STATE" "$STATE_GRID" "$RUN" "$GRID" \
  > "$RUN/regrid.log" 2>&1 || { echo "error: see $RUN/regrid.log" >&2; cat "$RUN/regrid.log" >&2; exit 1; }
tail -1 "$RUN/regrid.log"

echo "== 3. run (log: $RUN/run.log)"
start=$(date +%s)
set +e
timeout "$TIMEOUT" "$MPIRUN" -np "$NP" "$HERMES" -d "$RUN" restart mesh:file="$GRID" "$@" > "$RUN/run.log" 2>&1
rc=$?
set -e
echo "   exit code $rc after $(( $(date +%s) - start )) s"
[ $rc -eq 0 ] || { echo "error: the run did not finish; see $RUN/run.log and $RUN/BOUT.log.0" >&2; tail -5 "$RUN/run.log" >&2; exit $rc; }

echo "== 4. summary -> $RUN/summary.json"
"$PYTHON" "$HERE/summary.py" "$RUN" "$GRID" --json "$RUN/summary.json" > "$RUN/summary.txt"
"$PYTHON" - "$RUN/summary.json" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
c = d.get("change_per_ms", {})
print("   t = %.1f ms; change per ms over the last %.1f ms: electrons %.1e, neutrals %.1e, energy %.1e"
      % (1e3 * d["time_s"], c.get("interval_ms", 0), c.get("electron_content", 0),
         c.get("neutral_content", 0), c.get("stored_energy", 0)))
s = d["outer_midplane"]["separatrix"]
print("   outboard midplane separatrix: n_e %.3g m^-3, T_e %.0f eV, T_i %.0f eV" % (s["n_e"], s["T_e"], s["T_i"]))
p = d["power_W"]
print("   power [MW]: heating %.2f, sheath %.2f, impurity radiation %.2f"
      % (p["heating"] / 1e6, (p["sheath_electrons"] + p["sheath_ions"]) / 1e6, sum(p["impurity_radiation"].values()) / 1e6))
EOF
