"""Solver-recipe tests (us_runs/t23_recipe): cost and solution of every run.

usage (from us_runs/): python numerics/recipe_analyse.py out.json <reference run> <run> ...

For every run (see timing_analyse.py for the distance measure):
  time to 1 % / 0.1 %   wall time after which the state stays within 1 % (0.1 %)
                        of the run's own final state (only for runs that reached
                        the end time; otherwise None)
  steps, failures, Newton iterations per step
  diff_ref              difference of the final state from the reference run's final
                        state, per field and largest (volume-weighted 2-norms):
                        how much the change alters the solution
  t_end_ms, dist_ref_end  for runs stopped at the time limit: plasma time reached and
                        the distance of the last state from the reference final state
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import timing_analyse as ta  # noqa: E402


def series(run):
    """plasma time [ms], fields, volume weights (structured or flat run)"""
    import netCDF4
    import re
    inp = Path(run, "BOUT.inp").read_text()
    grid = re.search(r'^file\s*=\s*"([^"]+)"', inp, re.M).group(1)
    flat = "us_ncells" in netCDF4.Dataset(grid).variables
    t, wall, data, w, n = ta.flat_series(run, grid) if flat else ta.struct_series(run)
    return 1e3 * t, data, w


def main():
    out_path, ref_run, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
    _, ref_data, w = series(ref_run)
    ref_final = {f: a[-1] for f, a in ref_data.items()}
    t_end_nominal = None
    rows = []
    for run in runs:
        if not Path(run, "BOUT.dmp.0.nc").exists():
            continue
        res = ta.analyse(run)
        t_ms, data, wr = series(run)
        if t_end_nominal is None:
            t_end_nominal = t_ms[-1] if run == ref_run else None
        # reached the end time it was given (300 ms, or less for a test that is
        # steady long before): not stopped at the time limit, not aborted
        complete = "Run time" in Path(run, "run.log").read_text(errors="replace")
        diff = {}
        for f, a in data.items():
            if a.shape[1] != len(w):
                diff = None
                break
            diff[f] = float(np.sqrt(np.sum(w * (a[-1] - ref_final[f])**2) / np.sum(w * ref_final[f]**2)))
        w1, t1 = ta.time_to(res, 1e-2) if complete else (None, None)
        w3, t3 = ta.time_to(res, 1e-3) if complete else (None, None)
        log = Path(run, "run.log").read_text(errors="replace")
        # plasma time of the last step in the log (also for runs that stopped early),
        # and exit code and wall time from the queue log
        import re
        times = re.findall(r"Time: ([0-9.e+-]+), timestep", Path(run, "BOUT.log.0").read_text(errors="replace"))
        t_last_ms = 1e3 * float(times[-1]) / 9.5788e7 if times else 0.0
        qlog = Path(run).parent.parent / "rqueue.log"
        m = re.findall(r"end   %s exit=(\d+) wall=([0-9.]+)" % re.escape(run), qlog.read_text()) if qlog.exists() else []
        exit_code, wall_run = (int(m[-1][0]), float(m[-1][1])) if m else (None, None)
        reason = re.findall(r"Reason: ([A-Z_]+)", log)
        row = {"run": run, "name": Path(run).name, "complete": bool(complete),
               "exit_ok": "Run time" in log, "exit_code": exit_code, "wall_run_s": wall_run,
               "t_last_ms": t_last_ms, "abort_reason": reason[-1] if reason else None,
               "cells": res["cells"], "ranks": res["ranks"],
               "steps": res["steps"], "failures": res["failures"],
               "newton_per_step": res["newton_per_step"], "krylov_per_newton": res["krylov_per_newton"],
               "wall_total_s": res["wall_s"][-1], "t_end_ms": float(t_ms[-1]),
               "wall_1pct_s": w1, "t_1pct_ms": t1, "wall_01pct_s": w3, "t_01pct_ms": t3,
               "diff_ref": diff, "diff_ref_max": max(diff.values()) if diff else None,
               "t_ms": res["t_ms"], "wall_s": res["wall_s"], "distance": res["distance"]}
        if diff is not None:
            row["dist_ref_series"] = [float(max(np.sqrt(np.sum(w * (data[f][k] - ref_final[f])**2)
                                                         / np.sum(w * ref_final[f]**2)) for f in data))
                                      for k in range(len(t_ms))]
        rows.append(row)
        print("%-20s %s steps %5d fail %3d Newton/step %.2f Krylov/Newton %5.2f total %6.0f s | "
              "1%%: %s | 0.1%%: %s | t_end %.0f ms | final vs reference %s"
              % (row["name"], "done " if complete else "STOP ", row["steps"], row["failures"],
                 row["newton_per_step"] or 0, row["krylov_per_newton"] or 0, row["wall_total_s"],
                 "%5.0f s (%3.0f ms)" % (w1, t1) if w1 is not None else "   --        ",
                 "%5.0f s (%3.0f ms)" % (w3, t3) if w3 is not None else "   --        ",
                 row["t_end_ms"], "%.2e" % row["diff_ref_max"] if diff else "--"), flush=True)
    json.dump(rows, open(out_path, "w"))


if __name__ == "__main__":
    main()
