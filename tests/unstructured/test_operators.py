"""Tests of the flat-cell grid tools and the reference operators (no Hermes-3 run needed).

Run with:  python -m pytest tests/unstructured
      or:  python tests/unstructured/test_operators.py     (no pytest needed)
"""
import sys
from pathlib import Path

import numpy as np

try:
    import pytest
except ImportError:      # minimal stand-in so that the file runs on its own
    class pytest:        # noqa: N801
        class mark:      # noqa: N801
            @staticmethod
            def parametrize(names, values):
                def wrap(fn):
                    fn.cases = [v if isinstance(v, tuple) else (v,) for v in values]
                    return fn
                return wrap

        @staticmethod
        def approx(value, rel=1e-6):
            class Approx:
                def __eq__(self, other):
                    return abs(other - value) <= rel * abs(value)
            return Approx()

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "unstructured"))
import annulus  # noqa: E402
import flatgrid  # noqa: E402
import usops  # noqa: E402


def geometry(g):
    """Exact unit vectors of the annulus at the face mid-points"""
    xm = np.stack([g["face_Rmid"], g["face_Zmid"]], 1)
    rad = xm - np.array([annulus.R0, 0.0])
    er = rad / np.hypot(*rad.T)[:, None]
    bhat = np.stack([-er[:, 1], er[:, 0]], 1)                 # along the poloidal field
    Bp = annulus.PSI0 * np.hypot(*rad.T) / xm[:, 0]
    Bt = annulus.RBT / xm[:, 0]
    nrm = np.stack([g["face_nR"], g["face_nZ"]], 1)
    n_r = np.einsum("ij,ij->i", nrm, er)[:, None] * er        # radial part of the normal
    bn = np.einsum("ij,ij->i", nrm, bhat) * Bp / np.hypot(Bp, Bt)
    n_iso = nrm - bn[:, None] * bhat * (Bp / np.hypot(Bp, Bt))[:, None]
    return xm, n_r, n_iso


def cell_positions(g):
    return np.stack([g["R"], g["Z"]], 1)


@pytest.mark.parametrize("skew,jitter", [(0.0, 0.0), (0.4, 0.0), (0.4, 0.3)])
def test_grid_is_consistent(skew, jitter):
    g = annulus.build(8, 32, skew, jitter)
    out = flatgrid.check(g, verbose=False)
    assert out["max |sum of face fluxes| / sum |flux| per real cell"] < 1e-12
    assert out["max |sum of face normals * length| / perimeter per real cell"] < 1e-12
    assert out["closed flux tubes"] == 8
    assert (g["vertex_kind"] < 10).all()        # every vertex has a linear fit


@pytest.mark.parametrize("isotropic", [False, True])
def test_linear_field_is_exact_with_correction(isotropic):
    g = annulus.build(8, 32, skew=0.4, jitter=0.3)
    xm, n_r, n_iso = geometry(g)
    x = cell_positions(g)
    interior = g["face_kind"] == flatgrid.KIND_INTERIOR
    # faces whose vertices are both interior: the core boundary polygon is not straight
    clean = interior & (g["vertex_kind"][g["face_v1"]] != 1) & (g["vertex_kind"][g["face_v2"]] != 1)
    worst_two_point = 0.0
    for grad in (np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.6, -0.8])):
        f = x @ grad
        # with the face-mean geometry of the grid file the flow of a linear field is
        # (transmissibility) * (f_R - f_L) evaluated with the same normal: compare both
        # schemes with the two-point flow on the line through the face mid-point
        nvec = n_iso if isotropic else n_r
        exact = g["face_area"] * (nvec @ grad)
        two = usops.face_gradient(g, f, isotropic)
        cor = usops.face_gradient(g, f, isotropic, skew=True, limit=None)
        scale = g["face_area"]
        assert np.abs(cor - exact)[clean].max() / scale[clean].max() < 2e-2
        worst_two_point = max(worst_two_point, (np.abs(two - exact) / scale)[clean].max())
        assert np.percentile((np.abs(cor - exact) / scale)[clean], 90) < 5e-3
    assert worst_two_point > 0.2      # the two-point flow is wrong by tens of percent here


def face_flow_error(n, skew, jitter, isotropic=True):
    g = annulus.build(n, 4 * n, skew, jitter)
    xm, n_r, n_iso = geometry(g)
    x = cell_positions(g)
    k = 9.0
    f = np.sin(k * x[:, 0]) * np.cos(k * x[:, 1])
    grad = np.stack([k * np.cos(k * xm[:, 0]) * np.cos(k * xm[:, 1]),
                     -k * np.sin(k * xm[:, 0]) * np.sin(k * xm[:, 1])], 1)
    exact = g["face_area"] * np.einsum("ij,ij->i", n_iso if isotropic else n_r, grad)
    sel = (g["face_kind"] == flatgrid.KIND_INTERIOR) & (g["vertex_kind"][g["face_v1"]] == 0) & (g["vertex_kind"][g["face_v2"]] == 0)
    err = lambda a: np.sqrt(np.mean(((a - exact) / g["face_area"])[sel]**2)) / k
    return err(usops.face_gradient(g, f, isotropic)), err(usops.face_gradient(g, f, isotropic, skew=True, limit=None))


def test_corrected_flow_converges_on_a_skewed_grid():
    two, cor = zip(*(face_flow_error(n, skew=0.4, jitter=0.0) for n in (8, 16, 32)))
    order = np.log2(cor[0] / cor[2]) / 2
    assert order > 1.5, (cor, order)
    assert two[2] > 0.5 * two[0]          # the two-point flow does not converge
    assert cor[2] < 0.05 * two[2]


def test_two_point_flow_converges_without_skew():
    two, cor = zip(*(face_flow_error(n, skew=0.0, jitter=0.0) for n in (8, 16, 32)))
    assert np.log2(two[0] / two[2]) / 2 > 1.5
    assert np.log2(cor[0] / cor[2]) / 2 > 1.5


def test_operator_conserves():
    g = annulus.build(8, 32, skew=0.4, jitter=0.3)
    rng = np.random.default_rng(3)
    f = rng.uniform(1, 2, len(g["vol"]))
    a = rng.uniform(1, 2, len(g["vol"]))
    for op in (lambda: usops.div_a_grad_perp(g, a, f, skew=True),
               lambda: usops.div_a_grad_perp_upwind(g, a, f, isotropic=True, skew=True)):
        total = (op() * g["vol"])[:g["nreal"]].sum()
        # only core faces exchange with the outside: wall faces carry no flow
        core = g["face_kind"] == flatgrid.KIND_CORE
        grad = usops.face_gradient(g, f, True, skew=True)
        assert (grad[g["face_kind"] == flatgrid.KIND_WALL] == 0).all()
        assert abs(total) < np.abs(grad[core] * 2).sum() + 1e-12
    # interior faces cancel exactly: with the core faces removed nothing is left
    gi = dict(g)
    gi["face_tiso"] = np.where(g["face_kind"] == 0, g["face_tiso"], 0.0)
    gi["face_viso"] = np.where(g["face_kind"] == 0, g["face_viso"], 0.0)
    total = (usops.div_a_grad_perp(gi, a, f, skew=True) * g["vol"])[:g["nreal"]].sum()
    assert abs(total) < 1e-12


def test_limiter_keeps_the_direction_of_the_two_point_flow():
    g = annulus.build(8, 32, skew=0.6, jitter=0.4)
    rng = np.random.default_rng(5)
    f = rng.uniform(0, 1, len(g["vol"]))
    two = usops.face_gradient(g, f, True)
    lim = usops.face_gradient(g, f, True, skew=True, limit=1.0)
    assert (np.sign(lim) == np.sign(two)).all() or (lim[np.sign(lim) != np.sign(two)] == 0).all()
    assert (np.abs(lim) <= 2 * np.abs(two) + 1e-14).all()


def test_wall_loss_total():
    g = annulus.build(8, 32, skew=0.4, jitter=0.0)
    n = np.full(len(g["vol"]), 2.0)
    D = np.full(len(g["vol"]), 0.5)
    lam = 0.03
    loss = -(usops.wall_loss(g, D, n, 1 / lam) * g["vol"])[:g["nreal"]].sum()
    wall_area = g["face_area"][g["face_kind"] == flatgrid.KIND_WALL].sum()   # per radian
    assert loss == pytest.approx(0.5 * 2.0 * wall_area / lam, rel=1e-12)
    # the wall is a flux surface here, so the whole face area takes part
    exact_area = annulus.R0 * 2 * np.pi * annulus.R_OUT
    assert wall_area == pytest.approx(exact_area, rel=1e-2)


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        for case in getattr(fn, "cases", [()]):
            try:
                fn(*case)
                print("PASS  %s%s" % (name, case if case else ""))
            except AssertionError as e:
                failed += 1
                print("FAIL  %s%s  %s" % (name, case if case else "", e))
    sys.exit(1 if failed else 0)
