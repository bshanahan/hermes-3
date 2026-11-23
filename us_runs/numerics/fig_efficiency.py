"""Figure: cost of TCV-X21 runs to steady state.

usage (from us_runs/): python numerics/fig_efficiency.py <timing json> <out stem>
(the json is written by numerics/timing_analyse.py for the runs named below)

(a) distance to the final state against wall-clock time, from a cold start
(b) time step against plasma time for the same runs
(c) parallel speed-up for 50 ms continued from steady state, with the parts
    measured by PETSc's profiler
(d) wall time to steady state (within 1 %) against the number of cells
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

timing = {Path(r["run"]).name: r for r in json.load(open(sys.argv[1]))}

RUNS = [("struct", ps.SECONDARY, "-", r"orthogonal grid, 13 processes"),
        ("flat_two_point", ps.ORANGE, "-", r"GOAT grid, two-point flux"),
        ("flat_skew_plasma", ps.AQUA, "--", r"GOAT grid, correction for the plasma"),
        ("flat_default", ps.BLUE, "-", r"GOAT grid, correction for plasma and atoms (as run)")]


def first_below(r, level):
    d = np.asarray(r["distance"])
    above = np.flatnonzero(d >= level)
    k = above[-1] + 1 if len(above) else 0
    return (r["wall_s"][k], r["t_ms"][k]) if k < len(d) else (np.nan, np.nan)


def steps(run):
    txt = Path("t22_timing", run, "BOUT.log.0").read_text(errors="replace")
    m = re.findall(r"Time: ([0-9.e+-]+), timestep: ([0-9.e+-]+), nl iter: \d+, lin iter: \d+, reason: (-?\d+)", txt)
    t = np.array([float(a) for a, _, _ in m]) / OMEGA
    dt = np.array([float(b) for _, b, _ in m]) / OMEGA
    ok = np.array([int(c) > 0 for _, _, c in m])
    return t, dt, ok


def petsc_times(path):
    """time and count of each PETSc event of interest, and the Newton iterations of the run"""
    out = {}
    txt = Path(path).read_text()
    for line in txt.splitlines():
        p = line.split()
        if p and p[0] in ("SNESSolve", "MatFDColorApply", "MatLUFactorNum", "KSPSolve", "SNESFunctionEval"):
            out[p[0]] = (float(p[3]), int(p[1]))
    log = (Path(path).parent / "BOUT.log.0").read_text(errors="replace")
    out["newton"] = sum(int(v) for v in re.findall(r"nl iter: (\d+), lin iter: \d+, reason: \d", log))
    return out


def per_unit(parts):
    """seconds per Newton iteration (whole solve), per Jacobian and per LU factorisation"""
    return {"solve": parts["SNESSolve"][0] / parts["newton"],
            "jac": parts["MatFDColorApply"][0] / parts["MatFDColorApply"][1],
            "lu": parts["MatLUFactorNum"][0] / parts["MatLUFactorNum"][1]}


fig = plt.figure(figsize=(ps.WIDTH, 5.9))
top = fig.add_gridspec(1, 3, left=0.085, right=0.985, top=0.855, bottom=0.555, wspace=0.62)
bot = fig.add_gridspec(1, 1, left=0.385, right=0.985, top=0.44, bottom=0.075)
axs = np.array([[fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1])], [fig.add_subplot(top[0, 2]), fig.add_subplot(bot[0, 0])]])

# (a) approach to steady state
ax = axs[0, 0]
for run, color, ls, label in RUNS:
    r = timing.get(run)
    if r is None:
        continue
    w = np.asarray(r["wall_s"])
    d = np.asarray(r["distance"])
    keep = (w > 0) & (d > 0)
    ax.plot(w[keep], d[keep], color=color, ls=ls, lw=1.4)
    w1, _ = first_below(r, 1e-2)
    ax.plot([w1], [1e-2], marker="o", color=color, ms=4.5, zorder=5)
ax.axhline(1e-2, color=ps.RULE, lw=0.7, zorder=0)
ax.text(2.3, 1.25e-2, r"1\,\%", fontsize=7.5, color=ps.MUTED, va="bottom")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2, 1000)
ax.set_ylim(1e-5, 3)
ax.set_xlabel(r"wall-clock time [s]")
ax.set_ylabel(r"distance to the steady state")
ax.set_title(r"(a) From a cold start", loc="left", fontsize=9, pad=4)

# (b) time step
ax = axs[0, 1]
for run, color, ls, label in RUNS:
    if not Path("t22_timing", run, "BOUT.log.0").exists():
        continue
    t, dt, ok = steps(run)
    t0 = t[0] - dt[0]
    ax.plot(1e3 * (t[ok] - t0), dt[ok], color=color, ls=ls, lw=1.1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-5, 300)
ax.set_xlabel(r"plasma time [ms]")
ax.set_ylabel(r"time step [s]")
ax.set_title(r"(b) Time step", loc="left", fontsize=9, pad=4)

# (c) strong scaling, 50 ms continued from steady state, cost per unit of work
ax = axs[1, 0]
nps, units = [], []
for n in (1, 2, 4, 8, 12, 16):
    runs = [Path("t22_timing", "scale_np%d%s" % (n, suf), "petsc_log.txt") for suf in ("", "_rep")]
    runs = [f for f in runs if f.exists()]
    if runs:
        us = [per_unit(petsc_times(f)) for f in runs]
        nps.append(n)
        units.append({k: np.mean([u[k] for u in us]) for k in us[0]})
if nps:
    nps = np.array(nps)
    for key, color, marker, label in (("solve", ps.BLUE, "o", r"whole solve"),
                                      ("jac", ps.ORANGE, "s", r"Jacobian"),
                                      ("lu", ps.AQUA, "^", r"LU (MUMPS)")):
        sp = np.array([units[0][key] / u[key] for u in units])
        ax.plot(nps, sp, color=color, marker=marker, ms=4.5, label=label)
        print("speed-up %-6s" % key, " ".join("%d:%.2f" % (n, v) for n, v in zip(nps, sp)))
    ax.plot([1, 16], [1, 16], color=ps.RULE, lw=0.8, ls=":", zorder=0)
    ax.text(10.5, 6.2, r"ideal", fontsize=7.5, color=ps.MUTED, rotation=37)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8, 16])
    ax.set_xticklabels(["1", "2", "4", "8", "16"])
    ax.set_yticks([1, 2, 4, 8, 16])
    ax.set_yticklabels(["1", "2", "4", "8", "16"])
    ax.set_xlim(0.9, 18)
    ax.set_ylim(0.9, 18)
    ax.legend(loc="upper left", fontsize=7, handlelength=1.6, borderaxespad=0.3)
ax.set_xlabel(r"processes")
ax.set_ylabel(r"speed-up")
ax.set_title(r"(c) Parallel speed-up", loc="left", fontsize=9, pad=4)

# (d) wall-clock time to steady state for every configuration
ax = axs[1, 1]
BARS = [("struct", r"orthogonal", ps.SECONDARY),
        ("flat_two_point", r"GOAT, two-point", ps.BLUE),
        ("flat_skew_plasma", r"GOAT, correction for plasma", ps.BLUE),
        ("flat_default_np4", r"GOAT, as run", ps.BLUE),
        ("flat_default", r"GOAT, as run", ps.BLUE),
        ("flat_default_rep", r"GOAT, as run, repeated", ps.BLUE),
        ("flat_default_np16", r"GOAT, as run", ps.BLUE),
        ("flat_default_mc", r"GOAT, as run, MC limiter", ps.BLUE),
        ("flat_coarse", r"GOAT coarse, as run", ps.BLUE),
        ("flat_fine", r"GOAT fine, as run", ps.BLUE),
        ("full", r"GOAT whole plasma, as run", ps.BLUE)]
rows = [(timing[k], lab, c) for k, lab, c in BARS if k in timing]
y = np.arange(len(rows))[::-1]
for yy, (r, lab, c) in zip(y, rows):
    w1, _ = first_below(r, 1e-2)
    ax.barh(yy, w1 / 60, height=0.62, color=c, lw=0)
    ax.text(w1 / 60 + 0.2, yy, r"%.1f min, %d steps" % (w1 / 60, r["steps"]), va="center", fontsize=7.5, color=ps.INK)
ax.set_yticks(y)
ax.set_yticklabels([r"%s: %d cells, %d proc." % (lab, r["cells"], r["ranks"]) for r, lab, _ in rows], fontsize=7.5)
ax.set_xlim(0, 23.5)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_xlabel(r"wall-clock time to within 1\,\% of the steady state [min]")
ax.tick_params(axis="y", left=False, right=False, length=0)
ax.set_title(r"(d) All runs, from a cold start", loc="left", fontsize=9, pad=4)

handles = [Line2D([], [], color=c, ls=ls, lw=1.4, label=l) for _, c, ls, l in RUNS]
fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.53, 0.995), columnspacing=1.6, fontsize=8)
ps.save(fig, sys.argv[2])
print("wrote", sys.argv[2])
