"""Figure: convergence of the cross-field flow on skewed annuli.

usage (from us_runs/numerics): python fig_convergence.py <out stem>
Reads op_convergence.json (face flows, numpy reference operators) and
mms/mms_errors.json (Hermes-3 runs with manufactured solutions).
"""
import json
import sys

import numpy as np
from matplotlib.lines import Line2D

import paperstyle as ps

ps.setup()
plt = ps.plt

ops = json.load(open("op_convergence.json"))
mms = json.load(open("mms/mms_errors.json"))

SCHEMES = [("two_point", ps.ORANGE, "two-point flux"),
           ("bounded", ps.AQUA, "with skew correction and bound (as run)"),
           ("corrected", ps.BLUE, "with skew correction, no bound")]
GRIDS = [("quadrilaterals", "quad", "o", "-", "quadrilaterals"),
         ("triangles", "tri", "^", "--", "triangles")]
MMS_KEY = {"two_point": "two_point", "bounded": "default", "corrected": "unbounded"}

fig, axs = plt.subplots(2, 2, figsize=(ps.WIDTH, 5.4), sharex=True)
for col, field in enumerate(("smooth", "steep")):
    ax = axs[0, col]
    for key, color, _ in SCHEMES:
        for gname, _, marker, ls, _ in GRIDS:
            rows = ops[gname]
            h = np.array([r["h"] for r in rows]) * 1e3
            e = np.array([r[field][key] for r in rows])
            ax.plot(h, e, color=color, ls=ls, marker=marker, ms=4.5 if marker == "o" else 5)
    ax = axs[1, col]
    for key, color, _ in SCHEMES:
        for _, gtag, marker, ls, _ in GRIDS:
            rows = sorted([r for r in mms if r["grid"] == gtag and r["solution"] == field
                           and r["setting"] == MMS_KEY[key]], key=lambda r: r["h"])
            if not rows:
                continue
            h = np.array([r["h"] for r in rows]) * 1e3
            e = np.array([r["l2"] for r in rows])
            ax.plot(h, e, color=color, ls=ls, marker=marker, ms=4.5 if marker == "o" else 5)

for ax in axs.flat:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3, 110)
axs[0, 0].set_ylim(5e-5, 2)
axs[0, 1].set_ylim(5e-5, 2)
axs[1, 0].set_ylim(5e-5, 0.4)
axs[1, 1].set_ylim(5e-5, 0.4)
for ax in axs[1]:
    ax.set_xlabel(r"mean cell size $h$ [mm]")
axs[0, 0].set_ylabel(r"error of the face flux")
axs[1, 0].set_ylabel(r"error of the solution $n$")

ps.slope_triangle(axs[0, 0], 6, 2.2e-4, 2)
ps.slope_triangle(axs[0, 1], 6, 7e-5, 2)
ps.slope_triangle(axs[1, 0], 6, 6e-5, 2)
ps.slope_triangle(axs[1, 1], 6, 3.5e-4, 2)
ps.slope_triangle(axs[1, 0], 30, 2.2e-4, 1)

titles = [r"(a) Face fluxes, smooth field", r"(b) Face fluxes, edge-like field",
          r"(c) Hermes-3, manufactured solution, smooth", r"(d) Hermes-3, manufactured solution, edge-like"]
for ax, t in zip(axs.flat, titles):
    ax.set_title(t, loc="left", fontsize=9, pad=4)
axs[1, 1].text(0.97, 0.06, r"no bound: no steady state on" "\n" r"coarse grids (density $<0$)",
               transform=axs[1, 1].transAxes, ha="right", va="bottom", fontsize=7.5, color=ps.SECONDARY)

fig.tight_layout(h_pad=1.2, w_pad=1.5, rect=(0, 0, 1, 0.92))
handles = [Line2D([], [], color=c, lw=1.6, label=l) for _, c, l in SCHEMES]
leg = fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0), columnspacing=1.8)
handles = [Line2D([], [], color=ps.SECONDARY, ls=ls, marker=m, ms=4.5, label=l) for _, _, m, ls, l in GRIDS]
fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.965), columnspacing=1.8)
ps.save(fig, sys.argv[1])
print("wrote", sys.argv[1] + ".pdf/.png")
