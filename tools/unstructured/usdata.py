"""Read Hermes-3 output from a flat-cell (unstructured) run, in SI units."""
import numpy as np
import netCDF4
import flatgrid


class Run:
    """Output of a flat-cell run. Cells are returned in logical order: real
    cells first (flux tube by flux tube), then ghost cells, whatever the
    storage order and the number of processors."""

    def __init__(self, path, gridfile):
        import glob
        files = sorted(glob.glob(path + "/BOUT.dmp.*.nc"), key=lambda f: int(f.split(".")[-2]))
        self.files = [netCDF4.Dataset(f) for f in files]
        self.d = self.files[0]
        self.g = flatgrid.read(gridfile)
        self.myg = int(self.d["MYG"][...])
        self.n = self.g["nreal"]
        self.ntot = self.n + self.g["nghost"]
        self.t = np.asarray(self.d["t_array"][:]) / float(self.d["Omega_ci"][...])

    def raw(self, name):
        """Variable `name` as [time, cell] in output units and logical order"""
        m = self.myg
        a = np.concatenate([np.asarray(f[name][:, 0, m:f[name].shape[2] - m, 0]) for f in self.files], axis=1)
        return a[:, self.g["storage_index"]]

    def __call__(self, name, ghosts=False):
        """Variable `name` as [time, cell] in SI units"""
        a = self.raw(name)
        v = self.d[name]
        if "conversion" in v.ncattrs():
            a = a * v.getncattr("conversion")
        return a if ghosts else a[:, :self.n]

    def has(self, name):
        return name in self.d.variables

    def volume_integral(self, a):
        """Integral over the real cells [.. m^3]"""
        return 2 * np.pi * (a[..., :self.n] * self.g["vol"][:self.n]).sum(-1)


def plot_cells(ax, g, values, **kw):
    """Colour the real cells of grid `g` by `values`"""
    from matplotlib.collections import PolyCollection
    polys = []
    for c in g["cell_vertices"]:
        c = c[c >= 0]
        polys.append(np.stack([g["vertex_R"][c], g["vertex_Z"][c]], 1))
    pc = PolyCollection(polys, array=np.asarray(values)[:len(polys)], **kw)
    ax.add_collection(pc)
    ax.set_xlim(g["vertex_R"].min() - 0.02, g["vertex_R"].max() + 0.02)
    ax.set_ylim(g["vertex_Z"].min() - 0.02, g["vertex_Z"].max() + 0.02)
    ax.set_aspect("equal")
    return pc
