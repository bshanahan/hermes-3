"""Errors of the MMS runs (see mms_setup.py).

usage (from us_runs/numerics/mms): python ../mms_analyse.py out.json

For every finished run <grid type>_<solution>_n<N>_<setting>: the RMS
(volume-weighted) and largest error of the final density at the cell
centroids, relative to the range of the exact solution, and the relative
change over the last output (steadiness).
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "unstructured"))
import flatgrid  # noqa: E402
import mms_setup  # noqa: E402
import usdata  # noqa: E402

rows = []
for d in sorted(Path(".").glob("*_n*_*/")):
    m = re.match(r"(quad|tri)_(smooth|steep)_n(\d+)_([a-z0-9_]+)$", d.name)
    if not m or not (d / "BOUT.dmp.0.nc").exists() or "Run time" not in (d / "run.log").read_text(errors="replace"):
        continue
    kind, name, n, cfg = m.group(1), m.group(2), int(m.group(3)), m.group(4)
    grid = "grids/%s_%s_n%d.grd.nc" % (kind, name, n)
    r = usdata.Run(str(d), grid)
    g = r.g
    nreal = g["nreal"]
    N = r("Ni")
    x = np.stack([g["R"][:nreal], g["Z"][:nreal]], 1)
    ex, _ = mms_setup.exact(name, x)
    w = g["vol"][:nreal]
    rng = ex.max() - ex.min()
    err = N[-1][:nreal] - ex
    rows.append({"grid": kind, "solution": name, "n": n, "setting": cfg, "cells": int(nreal),
                 "h": float(np.sqrt(np.mean(g["area"][:nreal]))),
                 "l2": float(np.sqrt(np.sum(w * err**2) / np.sum(w)) / rng),
                 "max": float(np.abs(err).max() / rng),
                 "change": float(np.abs(N[-1][:nreal] - N[-2][:nreal]).max() / rng)})
    print("%-5s %-6s n=%-3d %-9s cells %6d  L2 %.3e  max %.3e  last change %.1e"
          % (kind, name, n, cfg, nreal, rows[-1]["l2"], rows[-1]["max"], rows[-1]["change"]))
json.dump(rows, open(sys.argv[1], "w"), indent=1)
