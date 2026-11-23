"""Figure: the same TCV-X21 case on the orthogonal grid and on the GOAT grid, with the
same solver recipe and 8 processes.

usage (from us_runs/): python numerics/fig_grids.py <struct json> <goat json> <out stem> [fastest|recipe]
(json written by numerics/recipe_analyse.py)

(a) distance to the run's own steady state against wall-clock time
(b) time step against plasma time
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as ps  # noqa: E402

ps.setup()
plt = ps.plt
OMEGA = 9.5788e7

struct = {r["name"]: r for r in json.load(open(sys.argv[1]))}
goat = {r["name"]: r for r in json.load(open(sys.argv[2]))}

KEY = sys.argv[4] if len(sys.argv) > 4 else "fastest"     # which recipe: fastest or recipe
RUNS = [(struct, KEY, ps.SECONDARY, "-", r"orthogonal grid (%d cells)"),
        (goat, KEY, ps.BLUE, "-", r"GOAT grid (%d cells), skew correction as used"),
        (goat, KEY + "_rep", ps.BLUE, ":", r"GOAT as used, repeated"),
        (goat, KEY + "_plasma_only" if KEY == "fastest" else "plasma_only", ps.AQUA, "--",
         r"GOAT, correction for the plasma only"),
        (goat, KEY + "_two_point" if KEY == "fastest" else "two_point", ps.ORANGE, "-.",
         r"GOAT, two-point flux")]


def steps(run):
    txt = Path(run, "BOUT.log.0").read_text(errors="replace")
    m = re.findall(r"Time: ([0-9.e+-]+), timestep: ([0-9.e+-]+), nl iter: \d+, lin iter: \d+, reason: (-?\d+)", txt)
    t = np.array([float(a) for a, _, _ in m]) / OMEGA
    dt = np.array([float(b) for _, b, _ in m]) / OMEGA
    ok = np.array([int(c) > 0 for _, _, c in m])
    return t, dt, ok


fig, (axa, axb) = plt.subplots(1, 2, figsize=(ps.WIDTH, 2.7))
handles = []
for src, name, color, ls, label in RUNS:
    r = src.get(name)
    if r is None:
        continue
    w = np.asarray(r["wall_s"])
    d = np.asarray(r["distance"])
    keep = (w > 0) & (d > 0)
    axa.plot(w[keep], d[keep], color=color, ls=ls, lw=1.3)
    if r["wall_1pct_s"] is not None:
        axa.plot([r["wall_1pct_s"]], [1e-2], marker="o", color=color, ms=4.5, mec="white", mew=0.6, zorder=5)
    t, dt, ok = steps(r["run"])
    axb.plot(1e3 * (t[ok] - (t[0] - dt[0])), dt[ok], color=color, ls=ls, lw=1.1)
    if label:
        handles.append(Line2D([], [], color=color, ls=ls, lw=1.3, label=label % r["cells"] if "%d" in label else label))
axa.axhline(1e-2, color=ps.RULE, lw=0.7, zorder=0)
axa.text(2.3, 1.25e-2, r"1\,\%", fontsize=7.5, color=ps.MUTED, va="bottom")
axa.set_xscale("log")
axa.set_yscale("log")
axa.set_xlim(2, 600)
axa.set_ylim(1e-5, 3)
axa.set_xlabel(r"wall-clock time [s]")
axa.set_ylabel(r"distance to steady state")
axa.set_title(r"(a) Approach to steady state", loc="left", fontsize=9, pad=4)
axb.set_xscale("log")
axb.set_yscale("log")
axb.set_xlim(1e-5, 300)
axb.set_xlabel(r"plasma time [ms]")
axb.set_ylabel(r"time step [s]")
axb.set_title(r"(b) Time step", loc="left", fontsize=9, pad=4)
fig.tight_layout(w_pad=2.0, rect=(0, 0, 1, 0.86))
fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02), fontsize=7.5,
           columnspacing=1.4, handlelength=2.2)
ps.save(fig, sys.argv[3])
print("wrote", sys.argv[3])
