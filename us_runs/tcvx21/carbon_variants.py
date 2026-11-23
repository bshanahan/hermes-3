"""Compare carbon runs on the TCV-X21 GOAT grid: content, radiation and, where a run
evolves them, the carbon ion temperatures and the parallel velocity of the atoms.

usage (from us_runs/): python tcvx21/carbon_variants.py <label>=<run dir> [...]
"""
import sys
sys.path.insert(0, "../tools/unstructured"); sys.path.insert(0, "tcvx21")
import numpy as np, usdata, json
import compare as c

twopi = 2*np.pi; qe = c.qe


def adas(kind):
    d = json.load(open("json_database/%s96_c.json" % kind))
    lt, ln, co = np.array(d["log_temperature"]), np.array(d["log_density"]), np.array(d["log_coeff"])

    def f(level, te, n_e):
        x = np.clip(np.log10(np.maximum(te, 1e-30)), lt[0], lt[-1]); y = np.clip(np.log10(np.maximum(n_e, 1.0)), ln[0], ln[-1])
        i = np.clip(np.searchsorted(lt, x) - 1, 0, len(lt) - 2); j = np.clip(np.searchsorted(ln, y) - 1, 0, len(ln) - 2)
        a = (x - lt[i])/(lt[i + 1] - lt[i]); b = (y - ln[j])/(ln[j + 1] - ln[j])
        t = co[level]
        return 10**((1 - a)*(1 - b)*t[i, j] + a*(1 - b)*t[i + 1, j] + (1 - a)*b*t[i, j + 1] + a*b*t[i + 1, j + 1])
    return f


plt_, prb = adas("plt"), adas("prb")
rows = []
for arg in sys.argv[1:]:
    label, run = arg.split("=")
    r = usdata.Run(run, c.GRID); g = r.g; n = r.n; it = -1; vol = g["vol"][:n]
    integ = lambda a: twopi*(a*vol).sum()
    ions = [s for s in ["c+"] + ["c+%d" % k for k in range(2, 7)] if r.has("N" + s)]
    Nd = r("Nd+")[it]; Te = r("Te")[it]; Td = r("Td+")[it]
    N = {s: r("N" + s)[it] for s in ["c"] + ions}
    ne = Nd + sum((k + 1)*N[s] for k, s in enumerate(ions))
    nz = sum(N[s] for s in ions)
    # radiation: line radiation of each state plus recombination radiation onto it
    states = ["c"] + ions
    rad = sum(integ(ne*N[s]*plt_(k, Te, ne)) for k, s in enumerate(states[:-1])) + sum(integ(ne*N[s]*prb(k, Te, ne)) for k, s in enumerate(states[1:]))
    psin = c.psin_RZ(g["R"][:n], g["Z"][:n])
    omp = np.flatnonzero((np.abs(g["Z"][:n]) < 0.012) & (g["R"][:n] > c.q["rmagx"])); omp = omp[np.argsort(psin[omp])]
    x = c.rsep_of_psin(psin[omp]); sep = lambda a: float(np.interp(0.0, x, a[omp]))
    hist = [sum(twopi*(r("N" + s)[k]*vol).sum() for s in ions) for k in range(len(r.t))]
    k5 = min(5, len(hist) - 1)
    out = {"label": label, "t [ms]": 1e3*r.t[it], "carbon ion content": integ(nz), "change of it over the last %d outputs [%%]" % k5: 100*(hist[-1]/hist[-1 - k5] - 1),
           "carbon atom content": integ(N["c"]), "sputtered [1/s]": integ(r("Sc_sputter")[it]) if r.has("Sc_sputter") else np.nan,
           "carbon radiation [kW]": rad/1e3, "separatrix n_e [1e18]": sep(ne)/1e18, "separatrix Te [eV]": sep(Te), "separatrix Ti [eV]": sep(Td),
           "separatrix carbon concentration [%]": 100*sep(nz/ne), "separatrix Zeff": sep((Nd + sum((k + 1)**2*N[s] for k, s in enumerate(ions)))/ne)}
    for k, s in enumerate(ions):
        out["content " + s] = integ(N[s])
    rows.append(out)
    print("== %s (%s), t = %.1f ms" % (label, run, 1e3*r.t[it]))
    if r.has("Tc+"):
        print("   carbon ion temperature against the deuterium ion temperature, weighted by the density of each state:")
        for s in ions:
            T = r("T" + s)[it]; w = N[s]*vol; ratio = T/Td
            o = np.argsort(ratio); cw = np.cumsum(w[o])/w.sum()
            pct = lambda p: ratio[o][np.searchsorted(cw, p)]
            print("     %-4s mean T %.1f eV (T_D+ weighted the same way %.1f eV) | T/T_D+: 5%% %.2f, median %.2f, 95%% %.2f | fraction of the state where T differs from T_D+ by more than 20%%: %.1f%%" % (
                s, (w*T).sum()/w.sum(), (w*Td).sum()/w.sum(), pct(0.05), pct(0.5), pct(0.95), 100*w[np.abs(ratio - 1) > 0.2].sum()/w.sum()))
    if r.has("NVc") and np.abs(r.raw("NVc")[it]).max() > 0:
        V = r("Vc")[it] if r.has("Vc") else r("NVc")[it]/(12*1.67262192e-27*np.maximum(N["c"], 1e6))
        w = N["c"]*vol; vth = np.sqrt(qe*np.maximum(r("Tc")[it], 0.01)/(12*1.67262192e-27))
        print("   carbon atoms: density-weighted |V_par| %.0f m/s, thermal speed sqrt(T/m) %.0f m/s, largest |V_par| %.0f m/s; parallel flux |n V| integrated over the volume %.2e m^-2 s^-1 m^3" % (
            (w*np.abs(V)).sum()/w.sum(), (w*vth).sum()/w.sum(), np.abs(V[N["c"] > 1e-3*N["c"].max()]).max(), integ(N["c"]*np.abs(V))))
keys = [k for k in rows[0] if k != "label"]
print()
print("%-46s" % "" + "".join("%16s" % r_["label"] for r_ in rows))
for k in keys:
    print("%-46s" % k + "".join("%16.4g" % r_.get(k, np.nan) for r_ in rows))
