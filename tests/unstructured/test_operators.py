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
    # without the smoothing term the bound is |two-point term|: the flow never reverses
    g = annulus.build(8, 32, skew=0.6, jitter=0.4)
    rng = np.random.default_rng(5)
    f = rng.uniform(0, 1, len(g["vol"]))
    two = usops.face_gradient(g, f, True)
    lim = usops.face_gradient(g, f, True, skew=True, limit=1.0, smooth=0.0)
    assert (np.sign(lim) == np.sign(two)).all() or (lim[np.sign(lim) != np.sign(two)] == 0).all()
    assert (np.abs(lim) <= 2 * np.abs(two) + 1e-14).all()


def test_bound_keeps_empty_cells_filling():
    # default bound (harmonic mean): the flow into a cell that holds nothing never
    # reverses, however large the correction; the correction never exceeds the bound
    g = annulus.build(8, 32, skew=0.6, jitter=0.4)
    rng = np.random.default_rng(7)
    f = rng.uniform(0, 1, len(g["vol"])) * (rng.uniform(0, 1, len(g["vol"])) > 0.3)
    L, R = g["face_L"], g["face_R"]
    two = usops.face_gradient(g, f, True)
    lim = usops.face_gradient(g, f, True, skew=True)
    flows = g["face_tiso"] > 0                  # wall faces carry nothing
    empty_L = flows & (f[L] == 0) & (f[R] > 0)  # flow into L is positive
    empty_R = flows & (f[R] == 0) & (f[L] > 0)  # flow into R is -gradient
    assert empty_L.sum() > 10 and empty_R.sum() > 10
    assert (lim[empty_L] > 0).all() and (lim[empty_R] < 0).all()
    corr = lim - two
    unb = usops.face_gradient(g, f, True, skew=True, limit=None) - two
    tau = g["face_tiso"]
    a, b = np.abs(f[L]), np.abs(f[R])
    M = np.sqrt(two**2 + (tau * 2 * a * b / (a + b + 1e-300))**2)
    assert (np.abs(corr) <= M + 1e-14).all() and (np.abs(corr) <= np.abs(unb) + 1e-14).all()


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


# ---- drifts

def annulus_B2(R, Z):
    """|B|^2 of the annulus field"""
    return (annulus.PSI0**2 * ((R - annulus.R0)**2 + Z**2) + annulus.RBT**2) / R**2


def exact_face_integrals(g, grad_phi, npts=200):
    """Per face, by quadrature along the face and with the normal from L to R:
    the flux of curl(b/B) and the E x B volume flow of a potential with a
    constant gradient. Written out in components, independently of the grid file:
      curl(g e_phi) = (-dg/dZ, (1/R) d(R g)/dR)          g = Bt / B^2
      v_E           = (Bt / B^2) e_phi x grad(phi) = g * (dphi/dZ, -dphi/dR)"""
    s = (np.arange(npts) + 0.5) / npts
    x1 = np.stack([g["vertex_R"][g["face_v1"]], g["vertex_Z"][g["face_v1"]]], 1)
    x2 = np.stack([g["vertex_R"][g["face_v2"]], g["vertex_Z"][g["face_v2"]]], 1)
    R = x1[:, [0]] + s * (x2[:, [0]] - x1[:, [0]])
    Z = x1[:, [1]] + s * (x2[:, [1]] - x1[:, [1]])
    gfun = lambda r, z: (annulus.RBT / r) / annulus_B2(r, z)
    h = 1e-6
    dg_dZ = (gfun(R, Z + h) - gfun(R, Z - h)) / (2 * h)
    dRg_dR = ((R + h) * gfun(R + h, Z) - (R - h) * gfun(R - h, Z)) / (2 * h)
    nR, nZ, length = g["face_nR"][:, None], g["face_nZ"][:, None], g["face_length"][:, None]
    curl = ((-dg_dZ * nR + dRg_dR / R * nZ) * R).mean(1) * length[:, 0]
    exb = (gfun(R, Z) * (grad_phi[1] * nR - grad_phi[0] * nZ) * R).mean(1) * length[:, 0]
    return curl, exb


def test_curl_flux_matches_quadrature_and_sums_to_zero():
    g = annulus.build(8, 32, skew=0.4, jitter=0.3)
    curl, _ = exact_face_integrals(g, np.zeros(2))
    notwall = g["face_kind"] != flatgrid.KIND_WALL
    scale = np.abs(curl[notwall]).max()
    assert np.abs(g["face_curl"] - curl)[notwall].max() < 1e-5 * scale
    assert (g["face_curl"][~notwall] == 0).all()
    # a sign error or a missing factor would show: the fluxes are not small
    assert scale > 1e-3
    # cells without a wall face: the fluxes sum to zero (div curl = 0)
    total = np.zeros(len(g["vol"]))
    np.add.at(total, g["face_L"], g["face_curl"])
    np.add.at(total, g["face_R"], -g["face_curl"])
    haswall = np.zeros(len(g["vol"]), bool)
    haswall[g["face_L"][~notwall]] = True
    real = np.arange(len(g["vol"])) < g["nreal"]
    assert np.abs(total[real & ~haswall]).max() < 1e-14
    assert np.abs(total[real & haswall]).max() > 1e-4      # closed wall: these do not


def exb_flow_error(n):
    g = annulus.build(n, 4 * n, skew=0.4, jitter=0.0)
    grad = np.array([0.6, -0.8])
    phi = cell_positions(g) @ grad
    _, exact = exact_face_integrals(g, grad)
    sel = (g["face_kind"] == flatgrid.KIND_INTERIOR) & (g["vertex_kind"][g["face_v1"]] == 0) & (g["vertex_kind"][g["face_v2"]] == 0)
    flow = usops.flow_exb(g, phi)
    return np.sqrt(np.mean(((flow - exact) / g["face_area"])[sel]**2)), np.sqrt(np.mean((exact / g["face_area"])[sel]**2))


def test_exb_flow_of_a_linear_potential_converges():
    err, size = zip(*(exb_flow_error(n) for n in (8, 16, 32)))
    # the vertex values of a linear potential are exact; what is left is the
    # variation of R Bt / B^2 along the face, second order in the face length
    assert err[0] < 0.02 * size[0]
    assert np.log2(err[0] / err[2]) / 2 > 1.5, err


def test_drift_divergence_conserves():
    g = annulus.build(8, 32, skew=0.4, jitter=0.3)
    rng = np.random.default_rng(7)
    f = rng.uniform(1, 2, len(g["vol"]))
    phi = rng.uniform(-1, 1, len(g["vol"]))
    # magnetic drift with one temperature: nothing crosses the wall, and the core
    # boundary carries its mean value, so the volume integral vanishes
    total = (usops.div_f_flow(g, f, 1.7 * g["face_curl"]) * g["vol"])[:g["nreal"]].sum()
    assert abs(total) < 1e-13
    # E x B: interior faces exchange nothing with the outside
    flow = usops.flow_exb(g, phi)
    wall = g["face_kind"] == flatgrid.KIND_WALL
    assert (usops.flow_exb(g, phi, through_wall=False)[wall] == 0).all()
    interior_flow = np.where(g["face_kind"] == 0, flow, 0.0)
    total = (usops.div_f_flow(g, f, interior_flow) * g["vol"])[:g["nreal"]].sum()
    assert abs(total) < 1e-13
    # into the wall only: what the plasma loses there is the flow times the cell value
    assert (flow[wall] >= 0).all() and (flow[wall] > 0).any() and (flow[wall] == 0).any()
    wall_flow = np.where(wall, flow, 0.0)
    lost = (usops.div_f_flow(g, f, wall_flow) * g["vol"])[:g["nreal"]].sum()
    assert lost == pytest.approx((flow * f[g["face_L"]])[wall].sum(), rel=1e-12)


def ring_crossing_drift_flow(n, second_order):
    """Net magnetic-drift flow of a radial profile p(r) across each boundary
    between rings of the annulus, relative to the flow through one face.
    Physically zero: p is constant on the boundary and the curl flux of a
    closed surface vanishes."""
    g = annulus.build(n, 4 * n, skew=0.4, jitter=0.0)
    r = np.hypot(g["R"] - annulus.R0, g["Z"])
    p = np.exp(-(r - annulus.R_IN) / 0.05)          # steep, like an edge pedestal
    flow = 1.7 * g["face_curl"]
    value = usops.drift_face_values(g, p, flow, second_order=second_order)
    L, R = g["face_L"], g["face_R"]
    tube = g["tube"]
    cross = (g["face_kind"] == 0) & (tube[L] != tube[R])
    out = np.where(tube[R] > tube[L], 1.0, -1.0)      # outward
    worst = 0.0
    for t in range(n - 1):
        m = cross & (np.minimum(tube[L], tube[R]) == t)
        net = (flow * value * out)[m].sum()
        scale = np.abs(flow * value)[m].mean()
        worst = max(worst, abs(net) / scale)
    return worst


def test_second_order_drift_flux_has_little_numerical_diffusion():
    up = [ring_crossing_drift_flow(n, False) for n in (8, 16, 32)]
    tvd = [ring_crossing_drift_flow(n, True) for n in (8, 16, 32)]
    # first-order upwinding moves a net amount across every ring, several
    # times the flow through one face, and refining does not reduce it
    # (measured 7.8, 7.9, 8.0): numerical cross-field diffusion
    assert min(up) > 5.0, up
    # the limited second-order values take most of it away (measured 0.63,
    # 0.22, 0.083), and the rest converges faster than first order
    assert tvd[0] < 0.1 * up[0] and tvd[2] < 0.02 * up[2], (up, tvd)
    assert np.log2(tvd[0] / tvd[2]) / 2 > 1.3, tvd


def test_second_order_face_value_of_a_linear_field():
    g = annulus.build(8, 32, skew=0.4, jitter=0.3)
    f = 3.0 + 0.7 * g["R"] - 0.4 * g["Z"]
    flow = np.ones(len(g["face_L"]))
    value = usops.drift_face_values(g, f, flow, second_order=True)
    L, R = g["face_L"], g["face_R"]
    inner = (g["face_kind"] == 0)
    mid = 0.5 * (f[L] + f[R])
    # between the two cell values, and close to their mean
    lo, hi = np.minimum(f[L], f[R]), np.maximum(f[L], f[R])
    assert ((value[inner] >= lo[inner] - 1e-12) & (value[inner] <= hi[inner] + 1e-12)).all()
    assert np.abs(value - mid)[inner].max() < 0.2 * np.abs(f[R] - f[L])[inner].max()


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
