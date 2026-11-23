"""Figure: numerical cross-field transport caused by parallel conduction.

usage (from us_runs/numerics): python fig_leakage.py <out stem>      (reads pollution.json)
"""
import json
import sys

import numpy as np
from matplotlib.lines import Line2D

import paperstyle as ps

ps.setup()
plt = ps.plt
data = json.load(open("pollution.json"))
k = np.array(data["ratios"])
use = k > 0

CASES = [("quadrilaterals, snapped", ps.BLUE, "-", "o", r"quadrilaterals, vertices on the flux surfaces"),
         ("quadrilaterals, psi error 3e-4", ps.ORANGE, "-", "s", r"quadrilaterals, vertex $\psi$ error $3\times10^{-4}$"),
         ("triangles, snapped", ps.AQUA, "--", "^", r"triangles, two-point parallel flux (as run)"),
         ("triangles, snapped, diamond parallel", ps.YELLOW, "-.", "v", r"triangles, parallel flux from the corrected gradient")]

fig, ax = plt.subplots(figsize=(4.6, 3.3))
ax.axvspan(1e7, 1e9, color="#e8e7e1", lw=0, zorder=0)
ax.text(1e8, 4.2, r"TCV-X21" "\n" r"edge electrons", ha="center", va="top", fontsize=7.5, color=ps.SECONDARY)
for key, color, ls, marker, label in CASES:
    for n, alpha, lw in ((16, 0.4, 1.0), (32, 1.0, 1.4)):
        e = np.array(data["cases"]["%s, n=%d" % (key.replace(", diamond parallel", ""), n)
                                   + (", diamond parallel" if "diamond" in key else "")])
        ax.plot(k[use], e[use], color=color, ls=ls, marker=marker if n == 32 else None, ms=4,
                alpha=alpha, lw=lw)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.7, 1.4e10)
ax.set_ylim(1e-3, 5)
ax.set_xlabel(r"$\kappa_\parallel / \kappa_\perp$")
ax.set_ylabel(r"error of the solution $T$")
handles = [Line2D([], [], color=c, ls=ls, marker=m, ms=4, label=l) for _, c, ls, m, l in CASES]
handles += [Line2D([], [], color=ps.SECONDARY, lw=1.0, alpha=0.4, label=r"lighter lines: half the resolution")]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.2), ncol=1, fontsize=7.5)
fig.tight_layout()
ps.save(fig, sys.argv[1])
print("wrote", sys.argv[1] + ".pdf/.png")
