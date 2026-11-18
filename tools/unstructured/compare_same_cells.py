"""Compare a structured Hermes-3 run with flat-cell runs on the same cells.

usage: python compare_same_cells.py <structured run dir> <flat run dir> <flat grid file> [time index]
"""
import sys
import numpy as np
import netCDF4
import usdata


def structured(path, name, it):
    d = netCDF4.Dataset(path + "/BOUT.dmp.0.nc")
    v = d[name]
    a = np.asarray(v[it, :, :, 0])
    if "conversion" in v.ncattrs():
        a = a * v.getncattr("conversion")
    myg = int(d["MYG"][...])
    if a.shape[1] > int(d["ny"][...]):
        a = a[:, myg:-myg]
    return a, np.asarray(d["t_array"][:]) / float(d["Omega_ci"][...])


def main(spath, fpath, gridfile, it=-1):
    r = usdata.Run(fpath, gridfile)
    f = netCDF4.Dataset(gridfile)
    n = r.n
    xi = np.asarray(f["us2d_structured_x"][0, :n]).astype(int)
    yj = np.asarray(f["us2d_structured_y"][0, :n]).astype(int)
    _, ts = structured(spath, "Nd+", it)
    print("structured t = %.3f ms, flat t = %.3f ms" % (1e3 * ts[it], 1e3 * r.t[it]))
    qe = 1.602176634e-19
    out = {}
    for name, label, scale in (("Nd+", "density", 1.0), ("Td+", "ion temperature", 1.0), ("Te", "electron temperature", 1.0),
                               ("Vd+", "ion velocity", 1.0)):
        a, _ = structured(spath, name, it)
        s = a[xi, yj]
        u = r(name)[it]
        if name == "Vd+":
            ref = np.abs(s).max()
            diff = np.abs(u - s) / ref
            txt = "|diff| / max|V|"
        else:
            diff = np.abs(u / s - 1)
            txt = "|ratio - 1|"
        k = np.argmax(diff)
        print("  %-22s %s: median %.4f  90%% %.4f  99%% %.4f  max %.4f at (x=%d, y=%d)" % (
            label, txt, np.median(diff), np.percentile(diff, 90), np.percentile(diff, 99), diff.max(), xi[k], yj[k]))
        out[name] = (s, u, diff)
    return out, xi, yj, r


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else -1)
