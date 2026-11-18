"""Record the headline numbers of a TCV-X21 GOAT-grid run, or check a run against a record.

usage (from us_runs/):
    python tcvx21/reference.py record <run dir> <file.json> ["note"]
    python tcvx21/reference.py check  <run dir> <file.json> [relative tolerance, default 0.02]

Recorded: midplane separatrix values, particle and power flows, and the
distance d of each profile to the forward- and reversed-field measurements
(as tcvx21/compare.py prints them). Needs the TCV-X21 dataset in tcvx21/data.
"""
import json
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np
import compare as c

PANELS = [("FHRP", "density"), ("FHRP", "electron_temp"), ("FHRP", "jsat"), ("TS", "density"), ("TS", "electron_temp"),
          ("LFS-LP", "density"), ("LFS-LP", "electron_temp"), ("LFS-LP", "jsat"), ("LFS-IR", "q_parallel"),
          ("HFS-LP", "density"), ("HFS-LP", "electron_temp"), ("HFS-LP", "jsat")]


def numbers(run):
    d = c.load(run); r, g, n = d["r"], d["g"], d["n"]
    twopi = 2*np.pi
    psin = d["psin"][:n]
    omp = np.flatnonzero((np.abs(g["Z"][:n]) < 0.012) & (g["R"][:n] > c.q["rmagx"])); omp = omp[np.argsort(psin[omp])]
    x = c.rsep_of_psin(psin[omp])
    sep = lambda a: float(np.interp(0.0, x, a[omp]))
    out = {"time [ms]": 1e3*float(d["t"]),
           "midplane separatrix n_e [m^-3]": sep(d["Ne"]), "midplane separatrix T_e [eV]": sep(d["Te"]),
           "midplane separatrix T_i [eV]": sep(d["Ti"]),
           "ion flow to walls through the sheath [1/s]": float(twopi*(d["gam"]*np.abs(d["fb"]))[d["wall"]].sum()),
           "sheath power [W]": float(twopi*(d["qpar"]*np.abs(d["fb"]))[d["wall"]].sum()),
           "ion content": float(r.volume_integral(d["Ne"])),
           "atom content": float(r.volume_integral(r("Nd", True)[-1]))}
    for diag, obs in PANELS:
        xs, vs = c.sim_profile(d, diag, obs)
        for field in ("forward", "reversed"):
            out["d %s %s %s" % (diag, obs, field)] = c.distance(xs, vs, *c.exp_profile(field, diag, obs))[0]
    return out


if __name__ == "__main__":
    mode, run, path = sys.argv[1:4]
    now = numbers(run)
    if mode == "record":
        json.dump({"note": sys.argv[4] if len(sys.argv) > 4 else "", "run": run, "values": now}, open(path, "w"), indent=1)
        print("recorded %d values in %s" % (len(now), path))
    else:
        tol = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
        ref = json.load(open(path))["values"]
        bad = 0
        for k, v in ref.items():
            rel = abs(now[k] - v)/max(abs(v), 1e-300)
            flag = "" if (rel <= tol or k.startswith("time")) else "   <-- differs"
            bad += bool(flag)
            print("%-46s recorded %12.5g  now %12.5g  (%+.2f%%)%s" % (k, v, now[k], 100*(now[k]/v - 1) if v else 0.0, flag))
        print("=> %s" % ("within %.0f%% of the record" % (100*tol) if not bad else "%d values differ by more than %.0f%%" % (bad, 100*tol)))
        sys.exit(1 if bad else 0)
