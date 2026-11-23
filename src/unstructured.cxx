#include "../include/unstructured.hxx"
#include "../include/hermes_build_config.hxx"

#include <bout/array.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

std::optional<Unstructured> Unstructured::instance;

namespace {
constexpr int halo_tag = 7741;

std::vector<int> getInts(Mesh* mesh, const std::string& name, int expected) {
  Array<int> a;
  if (mesh->get(a, name) != 0) {
    throw BoutException("Unstructured grid: '{}' missing from the mesh file", name);
  }
  if (a.size() != expected) {
    throw BoutException("Unstructured grid: '{}' has length {}, expected {}", name,
                        a.size(), expected);
  }
  return {a.begin(), a.end()};
}

std::vector<BoutReal> getReals(Mesh* mesh, const std::string& name, int expected,
                               BoutReal norm) {
  Array<BoutReal> a;
  if (mesh->get(a, name) != 0) {
    throw BoutException("Unstructured grid: '{}' missing from the mesh file", name);
  }
  if (a.size() != expected) {
    throw BoutException("Unstructured grid: '{}' has length {}, expected {}", name,
                        a.size(), expected);
  }
  std::vector<BoutReal> result(a.begin(), a.end());
  for (auto& v : result) {
    v /= norm;
  }
  return result;
}

/// Pick the entries of a whole-list array for a set of cells or faces
template <typename T>
std::vector<T> pick(const std::vector<T>& all, const std::vector<int>& which) {
  std::vector<T> result;
  result.reserve(which.size());
  for (const int i : which) {
    result.push_back(all[i]);
  }
  return result;
}

/// Value at a face seen from the upwind cell C, with a TVD limiter on an
/// unstructured grid (Darwish and Moukalled, Numer. Heat Transfer B 44, 2003):
/// the upwind-side difference a = 2 grad(f)_C . d_CD - (f_D - f_C) and the
/// downwind difference b = f_D - f_C are combined with van Albada's limiter,
/// f_face = f_C + s/2 with s = ab(a + b)/(a^2 + b^2) where ab > 0, else 0.
/// The switch at ab = 0 is smoothed so that Newton iterations see a
/// differentiable function. For a linear field s = b (the face mid-value);
/// s/b <= 1.21, so the face value stays between f_C and f_D.
inline BoutReal limitedFaceValue(BoutReal fc, BoutReal fd, BoutReal two_g_dot_d) {
  const BoutReal b = fd - fc;
  const BoutReal a = two_g_dot_d - b;
  const BoutReal scale = 1e-6 * (std::abs(fc) + std::abs(fd)) + 1e-300;
  const BoutReal ab = a * b;
  const BoutReal positive = 0.5 * (ab + std::sqrt(SQ(ab) + SQ(SQ(scale))));
  // (the last term keeps the denominator above zero when f_C = f_D = 0, where
  // SQ(scale) underflows)
  return fc + 0.5 * positive * (a + b) / (SQ(a) + SQ(b) + SQ(scale) + 1e-300);
}

/// Value of a reconstructed field on the side of a cell given by `side`
inline BoutReal faceValue(int side, BoutReal centre, BoutReal at_prev, BoutReal at_next) {
  return (side > 0) ? at_next : ((side < 0) ? at_prev : centre);
}
} // namespace

bool Unstructured::present(Mesh* mesh) { return mesh->sourceHasVar("us_ncells"); }

const Unstructured& Unstructured::get() {
  if (!instance) {
    throw BoutException("Unstructured grid used before it was initialised");
  }
  return *instance;
}

void Unstructured::initialise(Mesh* mesh, BoutReal rho_s0, BoutReal Bnorm) {
  int nproc = 1;
  int rank = 0;
  MPI_Comm_size(BoutComm::get(), &nproc);
  MPI_Comm_rank(BoutComm::get(), &rank);
  if (mesh->LocalNx != 1 or mesh->LocalNz != 1 or mesh->getNXPE() != 1) {
    throw BoutException(
        "Unstructured grids need nx = 1, nz = 1 and MXG = 0 (got nx = {}, nz = {})",
        mesh->LocalNx, mesh->LocalNz);
  }

  Unstructured us;
  int nreal = 0;
  int nghost = 0;
  int npad = 0;
  int nfaces = 0;
  mesh->get(nreal, "us_ncells");
  mesh->get(nghost, "us_nghost");
  mesh->get(npad, "us_npad");
  mesh->get(nfaces, "us_nfaces");
  const int ncells = nreal + nghost + npad; // The whole list

  us.nlocal = mesh->yend - mesh->ystart + 1;
  us.ys = mesh->ystart;
  us.global_start = rank * us.nlocal;
  us.length_norm = rho_s0;
  if (us.nlocal * nproc != ncells
      or mesh->getGlobalYIndexNoBoundaries(mesh->ystart) != us.global_start) {
    throw BoutException("Unstructured grid: {} cells cannot be split over {} processors "
                        "({} y points here, starting at {})",
                        ncells, nproc, us.nlocal,
                        mesh->getGlobalYIndexNoBoundaries(mesh->ystart));
  }

  //////////////////////////////////////////////////////////////////
  // The whole grid, as stored

  const BoutReal L = rho_s0;
  const auto all_vol = getReals(mesh, "us_vol", ncells, L * L * L);
  const auto all_B = getReals(mesh, "us_B", ncells, Bnorm);
  const auto all_area = getReals(mesh, "us_area", ncells, L * L);
  const auto all_bR = getReals(mesh, "us_bR", ncells, 1.0);
  const auto all_bZ = getReals(mesh, "us_bZ", ncells, 1.0);
  const auto all_cellR = getReals(mesh, "us_R", ncells, L);
  const auto all_cellZ = getReals(mesh, "us_Z", ncells, L);
  const auto all_prev = getInts(mesh, "us_prev", ncells);
  const auto all_next = getInts(mesh, "us_next", ncells);
  const auto all_kind = getInts(mesh, "us_kind", ncells);

  const auto all_faceL = getInts(mesh, "us_face_L", nfaces);
  const auto all_faceR = getInts(mesh, "us_face_R", nfaces);
  const auto all_face_kind = getInts(mesh, "us_face_kind", nfaces);
  const auto all_face_up = getInts(mesh, "us_face_up", nfaces);

  for (int f = 0; f < nfaces; f++) {
    const int cl = all_faceL[f];
    const int cr = all_faceR[f];
    const bool boundary = all_face_kind[f] != interior;
    if (cl < 0 or cl >= ncells or cr < 0 or cr >= ncells or all_kind[cl] != interior
        or (boundary != (all_kind[cr] != interior))
        or (boundary and (cl / us.nlocal != cr / us.nlocal))) {
      throw BoutException("Unstructured grid: face {} has inconsistent cells {}, {}", f,
                          cl, cr);
    }
  }

  // Ghost cells that close a flux tube: boundary faces crossed by the field,
  // on the open side of a cell at the end of its tube
  std::vector<int> all_ghost_prev(ncells, -1);
  std::vector<int> all_ghost_next(ncells, -1);
  for (int f = 0; f < nfaces; f++) {
    const int up = all_face_up[f];
    if (all_face_kind[f] == interior or up < 0) {
      continue;
    }
    const int c = all_faceL[f];
    if (up == all_prev[c] and all_next[c] < 0) {
      all_ghost_next[c] = all_faceR[f];
    } else if (up == all_next[c] and all_prev[c] < 0) {
      all_ghost_prev[c] = all_faceR[f];
    }
  }

  //////////////////////////////////////////////////////////////////
  // Cells each processor needs from the others: everything within two
  // faces of a cell it owns. That covers the face neighbours, their
  // neighbours along the flux tube and their ghost cells.

  // Skew correction: each face also needs the cells around its two vertices
  const bool has_stencils = mesh->sourceHasVar("us_face_sptr");
  us.skew_correction =
      Options::root()["unstructured"]["skew_correction"]
          .doc("Correct the cross-field operators for skewed cells, using values "
               "interpolated to the vertices? Needs a grid file written with the "
               "vertex stencils")
          .withDefault<bool>(has_stencils);
  us.skew_correction_mean =
      Options::root()["unstructured"]["skew_correction_neutrals"]
          .doc("Apply the skew correction to the diffusion of neutrals as well?")
          .withDefault<bool>(us.skew_correction);
  if (us.skew_correction_mean and !us.skew_correction) {
    throw BoutException("unstructured:skew_correction_neutrals needs skew_correction");
  }
  us.skew_limit =
      Options::root()["unstructured"]["skew_limit"]
          .doc("Largest size of the skew correction at a face, as a multiple of the "
               "two-point term (where the profile is under-resolved). Values up to 1 "
               "keep the flow into an empty cell from reversing")
          .withDefault<BoutReal>(1.0);
  us.skew_smooth =
      Options::root()["unstructured"]["skew_smooth"]
          .doc("The bound on the skew correction is at least this fraction of the "
               "field (mean of the two cells) times the transmissibility. Of order 1, "
               "it leaves the correction of resolved profiles untouched (second "
               "order); the first version used 0.05 with the arithmetic mean")
          .withDefault<BoutReal>(1.0);
  us.skew_harmonic =
      Options::root()["unstructured"]["skew_harmonic_mean"]
          .doc("Size of the field for skew_smooth by kind of field: harmonic mean of "
               "|f| for densities, temperatures and pressures (goes to zero when either "
               "cell empties, which keeps them positive), 1 for logarithms, root mean "
               "square for velocities and potentials? Otherwise the arithmetic mean of "
               "|f| for every field (first version; depends on the normalisation for "
               "logarithms)")
          .withDefault<bool>(true);
  us.skew_power =
      Options::root()["unstructured"]["skew_bound_power"]
          .doc("Sharpness p of the bound on the skew correction, C/(1+|C/M|^p)^(1/p)")
          .withDefault<BoutReal>(2.0);
  if (us.skew_power <= 0.0) {
    throw BoutException("unstructured:skew_bound_power must be positive");
  }
  if (us.skew_correction and !has_stencils) {
    throw BoutException("unstructured:skew_correction needs a grid file with vertex "
                        "stencils (us_face_sptr). Regenerate the grid");
  }
  std::vector<int> all_sptr(nfaces + 1, 0);
  std::vector<int> all_scell;
  std::vector<BoutReal> all_sw;
  us.has_drifts = has_stencils and mesh->sourceHasVar("us_face_drift");
  us.drift_second_order =
      Options::root()["unstructured"]["drift_second_order"]
          .doc("Drift fluxes with limited second-order face values? First-order "
               "upwind values add a numerical cross-field diffusion of about "
               "|v_drift| * (cell size) / 2")
          .withDefault<bool>(true);
  us.drift_second_order_closed_only =
      Options::root()["unstructured"]["drift_second_order_closed_only"]
          .doc("Use second-order drift face values only between cells of closed flux "
               "tubes (upwind values in the SOL and the divertor)? Most of the "
               "numerical diffusion that matters is across closed flux surfaces, and "
               "second-order values in the tenuous far SOL stall the implicit solver")
          .withDefault<bool>(true);
  // Normalised by Nnorm: read here in m^-3, divided below
  us.drift_second_order_density =
      Options::root()["unstructured"]["drift_second_order_density"]
          .doc("Density [m^-3] below which the drift fluxes fall back from second-order "
               "to upwind face values (weight n/(n + this)). Second-order values do not "
               "keep an implicit step positive, and in the tenuous far SOL they let "
               "cells drain. 0: second order everywhere")
          .withDefault<BoutReal>(1e17)
      / Options::root()["hermes"]["Nnorm"].withDefault<BoutReal>(1e19);
  us.exb_through_wall =
      Options::root()["unstructured"]["exb_through_wall"]
          .doc("Does the E x B drift carry plasma into the wall? If not, no drift "
               "flow crosses a wall face and the plasma it brings can only leave "
               "along the field. Needs a grid file written with vertex stencils "
               "for the wall faces")
          .withDefault<bool>(true);
  us.field_sign =
      Options::root()["unstructured"]["field_sign"]
          .doc("Sign of the toroidal field in the drifts: -1 reverses the field "
               "relative to the grid file. Parallel transport is not affected")
          .withDefault<BoutReal>(1.0);
  if (us.skew_correction or us.has_drifts) {
    all_sptr = getInts(mesh, "us_face_sptr", nfaces + 1);
    all_scell = getInts(mesh, "us_face_scell", all_sptr[nfaces]);
    all_sw = getReals(mesh, "us_face_sw", all_sptr[nfaces], 1.0);
  }

  std::vector<std::vector<int>> adjacent(ncells);
  for (int f = 0; f < nfaces; f++) {
    adjacent[all_faceL[f]].push_back(all_faceR[f]);
    adjacent[all_faceR[f]].push_back(all_faceL[f]);
    for (int k = all_sptr[f]; k < all_sptr[f + 1]; k++) {
      const int c = all_scell[k];
      if (c < 0 or c >= ncells) {
        throw BoutException("Unstructured grid: face {} has stencil cell {}", f, c);
      }
      adjacent[all_faceL[f]].push_back(c);
      adjacent[all_faceR[f]].push_back(c);
    }
  }
  const auto haloOf = [&](int proc) {
    const int lo = proc * us.nlocal;
    const int hi = lo + us.nlocal;
    std::vector<char> mark(ncells, 0);
    std::vector<int> first;
    for (int c = lo; c < hi; c++) {
      for (const int n : adjacent[c]) {
        if (mark[n] == 0) {
          mark[n] = 1;
          first.push_back(n);
        }
      }
    }
    for (const int c : first) {
      for (const int n : adjacent[c]) {
        mark[n] = 1;
      }
    }
    std::vector<int> halo; // Sorted, so grouped by owner
    for (int c = 0; c < ncells; c++) {
      if (mark[c] != 0 and (c < lo or c >= hi)) {
        halo.push_back(c);
      }
    }
    return halo;
  };

  const std::vector<int> halo = haloOf(rank);
  us.nall = us.nlocal + static_cast<int>(halo.size());
  std::vector<int> local(ncells, -1); // Local index of each cell in the list
  std::vector<int> cells(us.nall);    // And the reverse
  for (int c = 0; c < us.nlocal; c++) {
    local[us.global_start + c] = c;
    cells[c] = us.global_start + c;
  }
  for (std::size_t k = 0; k < halo.size(); k++) {
    local[halo[k]] = us.nlocal + static_cast<int>(k);
    cells[us.nlocal + k] = halo[k];
  }

  for (int proc = 0; proc < nproc; proc++) {
    if (proc == rank) {
      continue;
    }
    Exchange e{.rank = proc, .send = {}, .recv_start = 0, .recv_count = 0};
    for (const int c : haloOf(proc)) {
      if (c / us.nlocal == rank) {
        e.send.push_back(local[c]);
      }
    }
    for (std::size_t k = 0; k < halo.size(); k++) {
      if (halo[k] / us.nlocal == proc) {
        if (e.recv_count == 0) {
          e.recv_start = us.nlocal + static_cast<int>(k);
        }
        e.recv_count++;
      }
    }
    if (!e.send.empty() or e.recv_count > 0) {
      us.exchanges.push_back(e);
    }
  }

  //////////////////////////////////////////////////////////////////
  // Local cells

  const auto toLocal = [&](const std::vector<int>& all) {
    std::vector<int> result(us.nall);
    for (int c = 0; c < us.nall; c++) {
      const int other = all[cells[c]];
      result[c] = (other >= 0) ? local[other] : -1;
    }
    return result;
  };
  us.vol = pick(all_vol, cells);
  us.B = pick(all_B, cells);
  us.area = pick(all_area, cells);
  us.bR = pick(all_bR, cells);
  us.bZ = pick(all_bZ, cells);
  us.cellR = pick(all_cellR, cells);
  us.cellZ = pick(all_cellZ, cells);
  us.kind = pick(all_kind, cells);
  {
    // Closed flux tubes: every cell has a neighbour on both sides
    const auto all_tube = getInts(mesh, "us_tube", ncells);
    int ntubes = 0;
    for (int c = 0; c < ncells; c++) {
      ntubes = std::max(ntubes, all_tube[c] + 1);
    }
    std::vector<char> tube_closed(ntubes, 1);
    for (int c = 0; c < ncells; c++) {
      if (all_kind[c] == interior and (all_prev[c] < 0 or all_next[c] < 0)) {
        tube_closed[all_tube[c]] = 0;
      }
    }
    us.closed.assign(us.nall, 0);
    for (int c = 0; c < us.nall; c++) {
      const int t = all_tube[cells[c]];
      us.closed[c] = (all_kind[cells[c]] == interior and t >= 0) ? tube_closed[t] : 0;
    }
  }
  us.prev = toLocal(all_prev);
  us.next = toLocal(all_next);
  us.ghost_prev = toLocal(all_ghost_prev);
  us.ghost_next = toLocal(all_ghost_next);

  //////////////////////////////////////////////////////////////////
  // Local faces: those with a cell owned here

  std::vector<int> faces;
  for (int f = 0; f < nfaces; f++) {
    if (all_faceL[f] / us.nlocal == rank or all_faceR[f] / us.nlocal == rank) {
      faces.push_back(f);
    }
  }
  us.nfaces = static_cast<int>(faces.size());
  us.face_kind = pick(all_face_kind, faces);
  us.sideL = pick(getInts(mesh, "us_face_side_L", nfaces), faces);
  us.sideR = pick(getInts(mesh, "us_face_side_R", nfaces), faces);
  us.kpar = pick(getReals(mesh, "us_face_kpar", nfaces, L), faces);
  us.tperp = pick(getReals(mesh, "us_face_tperp", nfaces, L), faces);
  us.tiso = pick(getReals(mesh, "us_face_tiso", nfaces, L), faces);
  us.face_area = pick(getReals(mesh, "us_face_area", nfaces, L * L), faces);
  us.face_length = pick(getReals(mesh, "us_face_length", nfaces, L), faces);
  us.nR = pick(getReals(mesh, "us_face_nR", nfaces, 1.0), faces);
  us.nZ = pick(getReals(mesh, "us_face_nZ", nfaces, 1.0), faces);
  us.face_R = pick(getReals(mesh, "us_face_Rmid", nfaces, 1.0), faces);
  us.face_Z = pick(getReals(mesh, "us_face_Zmid", nfaces, 1.0), faces);
  // Flux [T m^2] / B [T]
  us.flux_over_B = pick(getReals(mesh, "us_face_flux", nfaces, L * L), faces);
  const auto faceB = pick(getReals(mesh, "us_face_B", nfaces, 1.0), faces);
  const auto faceBp = pick(getReals(mesh, "us_face_Bp", nfaces, 1.0), faces);
  if (us.has_drifts) {
    us.face_drift = pick(getReals(mesh, "us_face_drift", nfaces, L / Bnorm), faces);
    us.face_curl = pick(getReals(mesh, "us_face_curl", nfaces, L / Bnorm), faces);
    for (int i = 0; i < us.nfaces; i++) {
      us.face_drift[i] *= us.field_sign;
      us.face_curl[i] *= us.field_sign;
    }
  }
  if (us.skew_correction or us.has_drifts) {
    us.vperp = pick(getReals(mesh, "us_face_vperp", nfaces, L), faces);
    us.viso = pick(getReals(mesh, "us_face_viso", nfaces, L), faces);
    us.sptr.assign(1, 0);
    for (const int f : faces) {
      for (int k = all_sptr[f]; k < all_sptr[f + 1]; k++) {
        const int c = local[all_scell[k]];
        if (c < 0) {
          throw BoutException(
              "Unstructured grid: face {} has a stencil cell outside the halo", f);
        }
        us.scell.push_back(c);
        us.sw.push_back(all_sw[k]);
      }
      us.sptr.push_back(static_cast<int>(us.scell.size()));
    }
  }
  us.faceL.resize(us.nfaces);
  us.faceR.resize(us.nfaces);
  us.face_up.resize(us.nfaces);
  for (int i = 0; i < us.nfaces; i++) {
    const int f = faces[i];
    // Magnetic flux / (area * Bp) is the cosine between the face normal and the
    // poloidal field; the rest of the normal points across the flux surfaces
    const BoutReal bn_pol =
        std::min(1.0, std::abs(us.flux_over_B[i]) / (us.face_area[i] * std::max(faceBp[i], 1e-30)));
    us.wall_perp_area.push_back(
        (us.face_kind[i] == wall) ? us.face_area[i] * std::sqrt(1.0 - SQ(bn_pol)) : 0.0);
    us.flux_over_B[i] /= faceB[i];
    us.faceL[i] = local[all_faceL[f]];
    us.faceR[i] = local[all_faceR[f]];
    us.face_up[i] = (all_face_up[f] >= 0) ? local[all_face_up[f]] : -1;
    if (us.faceL[i] < 0 or us.faceR[i] < 0
        or (all_face_up[f] >= 0 and us.face_up[i] < 0)) {
      throw BoutException("Unstructured grid: face {} reaches a cell outside the halo",
                          f);
    }
  }

  // Cell volumes (per radian) and field strength into the BOUT++ coordinates,
  // so that J * dx * dy * dz is the cell volume
  auto* coord = mesh->getCoordinates();
  Field2D J{1.0, mesh};
  Field2D Bxy{1.0, mesh};
  J.allocate();
  Bxy.allocate();
  for (int c = 0; c < us.nlocal; c++) {
    J(0, us.ys + c) = us.vol[c];
    Bxy(0, us.ys + c) = us.B[c];
  }
  coord->setJ(J, false);
  coord->setBxy(Bxy, false);

  us.second_order =
      Options::root()["unstructured"]["second_order"]
          .doc("Reconstruct face values along flux tubes with a slope limiter? "
               "First order upwind if false")
          .withDefault<bool>(true);

  output.write("Unstructured grid: {} cells, {} ghost cells, {} faces in total; "
               "{} cells, {} halo cells and {} faces on this processor. "
               "Skew correction {} (bound: limit {}, smooth {}, {} scale, power {})\n",
               nreal, nghost, nfaces, us.nlocal, halo.size(), us.nfaces,
               us.skew_correction ? "on" : "off", us.skew_limit, us.skew_smooth,
               us.skew_harmonic ? "per-field" : "arithmetic", us.skew_power);
  instance = std::move(us);
}

Unstructured::Buffer Unstructured::gather(const Field3D& f) const {
  Buffer buf(nall);
  for (int c = 0; c < nlocal; c++) {
    buf[c] = f(0, ys + c, 0);
  }
  if (exchanges.empty()) {
    return buf;
  }
  std::vector<MPI_Request> requests;
  requests.reserve(2 * exchanges.size());
  std::vector<Buffer> outgoing(exchanges.size());
  for (const auto& e : exchanges) {
    if (e.recv_count > 0) {
      requests.emplace_back();
      MPI_Irecv(&buf[e.recv_start], e.recv_count, MPI_DOUBLE, e.rank, halo_tag,
                BoutComm::get(), &requests.back());
    }
  }
  for (std::size_t k = 0; k < exchanges.size(); k++) {
    const auto& e = exchanges[k];
    if (e.send.empty()) {
      continue;
    }
    outgoing[k].resize(e.send.size());
    for (std::size_t j = 0; j < e.send.size(); j++) {
      outgoing[k][j] = buf[e.send[j]];
    }
    requests.emplace_back();
    MPI_Isend(outgoing[k].data(), static_cast<int>(e.send.size()), MPI_DOUBLE, e.rank,
              halo_tag, BoutComm::get(), &requests.back());
  }
  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  return buf;
}

void Unstructured::reconstruct(const Buffer& f, Buffer& at_prev, Buffer& at_next) const {
  hermes::Limiter cellboundary;
  at_prev.resize(nall);
  at_next.resize(nall);
  for (int c = 0; c < nall; c++) {
    const BoutReal fc = f[c];
    // Neighbours along the tube. A ghost cell holds the value at its face, so
    // the equivalent guard-cell value is the reflection 2 * ghost - centre
    const bool has_m = prev[c] >= 0 or ghost_prev[c] >= 0;
    const bool has_p = next[c] >= 0 or ghost_next[c] >= 0;
    if (second_order and kind[c] == interior and has_m and has_p) {
      const BoutReal fm_ = (prev[c] >= 0) ? f[prev[c]] : 2. * f[ghost_prev[c]] - fc;
      const BoutReal fp_ = (next[c] >= 0) ? f[next[c]] : 2. * f[ghost_next[c]] - fc;
      FV::Stencil1D s{.c = fc, .m = fm_, .p = fp_};
      cellboundary(s);
      at_prev[c] = s.L;
      at_next[c] = s.R;
    } else {
      at_prev[c] = at_next[c] = fc;
    }
  }
}

Field3D Unstructured::Div_par_mod(const Field3D& f_in, const Field3D& v_in,
                                  const Field3D& wave_speed) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer f = gather(f_in);
  const Buffer v = gather(v_in);
  const Buffer wave = gather(wave_speed);
  Buffer fm, fp, vm, vp;
  reconstruct(f, fm, fp);
  reconstruct(v, vm, vp);

  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    if (face_kind[i] != interior) {
      // Ghost cell holds the boundary values
      result(0, ys + L, 0) += f[R] * v[R] * flux_over_B[i] / vol[L];
      continue;
    }
    const BoutReal fL = faceValue(sideL[i], f[L], fm[L], fp[L]);
    const BoutReal fR = faceValue(sideR[i], f[R], fm[R], fp[R]);
    const BoutReal vL = faceValue(sideL[i], v[L], vm[L], vp[L]);
    const BoutReal vR = faceValue(sideR[i], v[R], vm[R], vp[R]);
    const BoutReal amax = std::max({wave[L], wave[R], std::abs(v[L]), std::abs(v[R])});
    // Flux along b. The upstream state is L if b points from L to R
    const BoutReal flux = (flux_over_B[i] >= 0.0)
                              ? 0.5 * fL * (vL + amax) + 0.5 * fR * (vR - amax)
                              : 0.5 * fR * (vR + amax) + 0.5 * fL * (vL - amax);
    const BoutReal flow = flux * flux_over_B[i]; // From L to R
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_par_fvv(const Field3D& f_in, const Field3D& v_in,
                                  const Field3D& wave_speed, bool fixflux) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer f = gather(f_in);
  const Buffer v = gather(v_in);
  const Buffer wave = gather(wave_speed);
  Buffer fm, fp, vm, vp;
  reconstruct(f, fm, fp);
  reconstruct(v, vm, vp);

  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    if (face_kind[i] != interior) {
      const BoutReal nb = f[R]; // Boundary (mid-point) values
      const BoutReal vb = v[R];
      const bool at_next = ghost_next[L] == R;
      const bool at_prev = ghost_prev[L] == R;
      BoutReal flux; // Along b
      if (fixflux or not(at_next or at_prev)) {
        flux = nb * vb * vb;
      } else {
        // Cell-edge values, with a penalty that damps differences in velocity
        const BoutReal fe = at_next ? fp[L] : fm[L];
        const BoutReal ve = at_next ? vp[L] : vm[L];
        const BoutReal amax =
            std::max({wave[L], std::abs(v[L]), std::abs(2. * vb - v[L])});
        const BoutReal sign = (flux_over_B[i] >= 0.0) ? 1.0 : -1.0;
        flux = fe * ve * ve + sign * amax * nb * (ve - vb);
      }
      result(0, ys + L, 0) += flux * flux_over_B[i] / vol[L];
      continue;
    }
    const BoutReal fL = faceValue(sideL[i], f[L], fm[L], fp[L]);
    const BoutReal fR = faceValue(sideR[i], f[R], fm[R], fp[R]);
    const BoutReal vL = faceValue(sideL[i], v[L], vm[L], vp[L]);
    const BoutReal vR = faceValue(sideR[i], v[R], vm[R], vp[R]);
    const BoutReal amax = std::max({wave[L], wave[R], std::abs(v[L]), std::abs(v[R])});
    const BoutReal flux =
        (flux_over_B[i] >= 0.0)
            ? 0.5 * fL * (vL + amax) * vL + 0.5 * fR * (vR - amax) * vR
            : 0.5 * fR * (vR + amax) * vR + 0.5 * fL * (vL - amax) * vL;
    const BoutReal flow = flux * flux_over_B[i];
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_par(const Field3D& f_in) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer f = gather(f_in);
  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    if (face_kind[i] != interior) {
      result(0, ys + L, 0) += f[R] * flux_over_B[i] / vol[L];
      continue;
    }
    const BoutReal flow = 0.5 * (f[L] + f[R]) * flux_over_B[i];
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Grad_par(const Field3D& f_in) const {
  // b.Grad(f) = Div(b f) - f Div(b)
  Field3D result{zeroFrom(f_in)};
  const Buffer f = gather(f_in);
  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    if (face_kind[i] != interior) {
      result(0, ys + L, 0) += (f[R] - f[L]) * flux_over_B[i] / vol[L];
      continue;
    }
    const BoutReal half = 0.5 * (f[R] - f[L]) * flux_over_B[i];
    if (L < nlocal) {
      result(0, ys + L, 0) += half / vol[L];
    }
    if (R < nlocal) {
      result(0, ys + R, 0) += half / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_par_K_Grad_par(const Field3D& K_in, const Field3D& f_in) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer K = gather(K_in);
  const Buffer f = gather(f_in);
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] != interior) {
      continue;
    }
    const int L = faceL[i];
    const int R = faceR[i];
    const BoutReal flow = 0.5 * (K[L] + K[R]) * kpar[i] * (f[R] - f[L]); // Into L
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_a_Grad_perp_upwind(const Field3D& a_in, const Field3D& f_in,
                                             bool isotropic, Field3D* boundary_flow,
                                             Scale scale) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer a = gather(a_in);
  const Buffer f = gather(f_in);
  const auto& trans = isotropic ? tiso : tperp;
  for (int i = 0; i < nfaces; i++) {
    if (trans[i] == 0.0) {
      continue;
    }
    const int L = faceL[i];
    const int R = faceR[i];
    BoutReal gradient = trans[i] * (f[R] - f[L]);
    if (skew_correction) {
      gradient = correctedGradient(gradient, trans[i], isotropic ? viso[i] : vperp[i], f, i,
                                   scale);
    }
    // Flow is from high to low f: coefficient from the cell it leaves
    const BoutReal flow = gradient * ((gradient > 0) ? a[R] : a[L]); // Into L
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (boundary_flow != nullptr and face_kind[i] != interior) {
      (*boundary_flow)(0, ys + R, 0) -= flow; // Out of the plasma; both cells owned here
    }
    if (face_kind[i] == interior and R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_a_Grad_perp(const Field3D& a_in, const Field3D& f_in,
                                      Scale scale) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer a = gather(a_in);
  const Buffer f = gather(f_in);
  for (int i = 0; i < nfaces; i++) {
    if (tiso[i] == 0.0) {
      continue;
    }
    const int L = faceL[i];
    const int R = faceR[i];
    BoutReal gradient = tiso[i] * (f[R] - f[L]);
    if (skew_correction_mean) {
      gradient = correctedGradient(gradient, tiso[i], viso[i], f, i, scale);
    }
    const BoutReal flow = 0.5 * (a[L] + a[R]) * gradient; // Into L
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (face_kind[i] == interior and R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

std::vector<BoutReal> Unstructured::flowExB(const Field3D& phi_in) const {
  if (!has_drifts) {
    throw BoutException("Unstructured grid: the grid file has no drift geometry "
                        "(us_face_drift). Regenerate the grid");
  }
  const Buffer phi = gather(phi_in);
  std::vector<BoutReal> flow(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    flow[i] = face_drift[i] * vertexDifference(phi, i);
    if (face_kind[i] == wall) {
      // Plasma drifts into the wall but nothing drifts out of it. The cut-off
      // is smooth, with the same width as the upwind switch in Div_f_flow.
      const BoutReal width = 1e-3 * std::abs(face_drift[i]);
      flow[i] = exb_through_wall ? 0.5 * (flow[i] + std::sqrt(SQ(flow[i]) + SQ(width)))
                                 : 0.0;
    }
  }
  return flow;
}

BoutReal Unstructured::coreMean(const Buffer& f) const {
  // Core faces and their ghost cells are owned by one processor each
  BoutReal local[2] = {0.0, 0.0};
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] == core and faceL[i] < nlocal) {
      local[0] += face_area[i] * f[faceR[i]];
      local[1] += face_area[i];
    }
  }
  BoutReal global[2] = {0.0, 0.0};
  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, BoutComm::get());
  return (global[1] > 0.0) ? global[0] / global[1] : 0.0;
}

std::vector<BoutReal> Unstructured::flowMagneticDrift(const Field3D& T_in,
                                                      BoutReal charge) const {
  if (!has_drifts) {
    throw BoutException("Unstructured grid: the grid file has no drift geometry "
                        "(us_face_curl). Regenerate the grid");
  }
  const Buffer T = gather(T_in);
  const BoutReal Tcore = coreMean(T);
  std::vector<BoutReal> flow(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] == wall) {
      continue;
    }
    const BoutReal Tface =
        (face_kind[i] == interior) ? 0.5 * (T[faceL[i]] + T[faceR[i]]) : Tcore;
    flow[i] = Tface * face_curl[i] / charge;
  }
  return flow;
}

void Unstructured::cellGradient(const Field3D& like, const Buffer& f, Buffer& gR,
                                Buffer& gZ) const {
  Buffer sR(nall, 0.0), sZ(nall, 0.0);
  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    const bool boundary = face_kind[i] != interior;
    const BoutReal face_value = boundary ? f[R] : 0.5 * (f[L] + f[R]);
    sR[L] += face_value * nR[i] * face_length[i];
    sZ[L] += face_value * nZ[i] * face_length[i];
    if (!boundary) {
      sR[R] -= face_value * nR[i] * face_length[i];
      sZ[R] -= face_value * nZ[i] * face_length[i];
    }
  }
  Field3D fR{zeroFrom(like)}, fZ{zeroFrom(like)};
  for (int c = 0; c < nlocal; c++) {
    if (kind[c] == interior) {
      fR(0, ys + c, 0) = sR[c] / area[c];
      fZ(0, ys + c, 0) = sZ[c] / area[c];
    }
  }
  gR = gather(fR);
  gZ = gather(fZ);
}

std::vector<BoutReal> Unstructured::secondOrderWeight(const Field3D& n_in) const {
  std::vector<BoutReal> weight(nfaces, 1.0);
  if (drift_second_order_closed_only) {
    for (int i = 0; i < nfaces; i++) {
      if (closed[faceL[i]] == 0 or closed[faceR[i]] == 0) {
        weight[i] = 0.0;
      }
    }
  }
  if (drift_second_order_density <= 0.0) {
    return weight;
  }
  const Buffer n = gather(n_in);
  for (int i = 0; i < nfaces; i++) {
    if (weight[i] == 0.0) {
      continue;
    }
    const BoutReal nL = std::max(n[faceL[i]], 0.0);
    const BoutReal nR = std::max(n[faceR[i]], 0.0);
    const BoutReal nf = 2.0 * nL * nR / (nL + nR + 1e-300);
    weight[i] *= nf / (nf + drift_second_order_density);
  }
  return weight;
}

Field3D Unstructured::Div_f_flow(const Field3D& f_in, const std::vector<BoutReal>& flow,
                                 Field3D* boundary_flow,
                                 const std::vector<BoutReal>* weight) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer f = gather(f_in);
  const BoutReal fcore = coreMean(f);
  Buffer gR, gZ;
  if (drift_second_order) {
    cellGradient(f_in, f, gR, gZ);
  }
  for (int i = 0; i < nfaces; i++) {
    if (flow[i] == 0.0) {
      continue;
    }
    const int L = faceL[i];
    const int R = faceR[i];
    // From L to R. Every core face carries the mean of the boundary values,
    // for flow in either direction. The core boundary is a flux surface: with
    // one value on it the magnetic drifts take as much out as they bring in,
    // so that no net current or particle flow crosses it. (Upwinding, or the
    // local boundary values, let a net current through, which the potential
    // can only balance with an enormous radial electric field.)
    // Upwind value, switching side smoothly where the flow reverses: a sharp
    // switch stalls the Newton iterations of the implicit solver. The width of
    // the switch is the flow of a potential difference of 1e-3 (normalised)
    // between the vertices of the face.
    const BoutReal width = 1e-3 * std::abs(face_drift[i]);
    const BoutReal wL = 0.5 * (1.0 + flow[i] / std::sqrt(SQ(flow[i]) + SQ(width)));
    BoutReal value;
    if (face_kind[i] == core) {
      value = fcore;
    } else if (face_kind[i] == wall) {
      value = f[L]; // Flow through a wall face is outwards: the cell's value
    } else if (drift_second_order) {
      const BoutReal dR = cellR[R] - cellR[L];
      const BoutReal dZ = cellZ[R] - cellZ[L];
      const BoutReal fromL = limitedFaceValue(f[L], f[R], 2.0 * (gR[L] * dR + gZ[L] * dZ));
      const BoutReal fromR = limitedFaceValue(f[R], f[L], -2.0 * (gR[R] * dR + gZ[R] * dZ));
      const BoutReal upwind = wL * f[L] + (1.0 - wL) * f[R];
      const BoutReal second = wL * fromL + (1.0 - wL) * fromR;
      value = (weight == nullptr) ? second : upwind + (*weight)[i] * (second - upwind);
    } else {
      value = wL * f[L] + (1.0 - wL) * f[R];
    }
    const BoutReal out = flow[i] * value;
    if (face_kind[i] != interior and boundary_flow != nullptr) {
      (*boundary_flow)(0, ys + R, 0) += out; // A boundary face and its cells are owned here
    }
    if (L < nlocal) {
      result(0, ys + L, 0) += out / vol[L];
    }
    if (face_kind[i] == interior and R < nlocal) {
      result(0, ys + R, 0) -= out / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::wallLoss(const Field3D& a, const Field3D& f, BoutReal inv_length,
                               Field3D* flow) const {
  Field3D result{zeroFrom(f)};
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] != wall or wall_perp_area[i] == 0.0) {
      continue;
    }
    const int y = ys + faceL[i]; // A wall face and both its cells are owned here
    const BoutReal out = a(0, y, 0) * f(0, y, 0) * wall_perp_area[i] * inv_length;
    result(0, y, 0) -= out / vol[faceL[i]];
    if (flow != nullptr) {
      (*flow)(0, ys + faceR[i], 0) += out;
    }
  }
  return result;
}

Field3D Unstructured::gradPerpSq(const Field3D& f_in) const {
  const Buffer f = gather(f_in);
  // Green-Gauss gradient in the poloidal plane
  Buffer gR(nall, 0.0), gZ(nall, 0.0);
  for (int i = 0; i < nfaces; i++) {
    const int L = faceL[i];
    const int R = faceR[i];
    const bool boundary = face_kind[i] != interior;
    const BoutReal face_value = boundary ? f[R] : 0.5 * (f[L] + f[R]);
    const BoutReal wR = face_value * nR[i] * face_length[i];
    const BoutReal wZ = face_value * nZ[i] * face_length[i];
    gR[L] += wR;
    gZ[L] += wZ;
    if (!boundary) {
      gR[R] -= wR;
      gZ[R] -= wZ;
    }
  }
  Field3D result{zeroFrom(f_in)};
  for (int c = 0; c < nlocal; c++) {
    if (kind[c] != interior) {
      continue;
    }
    const BoutReal dR = gR[c] / area[c];
    const BoutReal dZ = gZ[c] / area[c];
    const BoutReal dpar = bR[c] * dR + bZ[c] * dZ;
    result(0, ys + c, 0) = std::max(dR * dR + dZ * dZ - dpar * dpar, 0.0);
  }
  // Ghost cells take the value of their neighbour
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] != interior) {
      result(0, ys + faceR[i], 0) = result(0, ys + faceL[i], 0);
    }
  }
  return result;
}

void Unstructured::setGhosts(Field3D& f, std::optional<BoutReal> core_value) const {
  for (int i = 0; i < nfaces; i++) {
    if (face_kind[i] == interior) {
      continue;
    }
    f(0, ys + faceR[i], 0) =
        (face_kind[i] == core and core_value) ? *core_value : f(0, ys + faceL[i], 0);
  }
}

void Unstructured::zeroGhosts(Field3D& f) const {
  for (int c = 0; c < nlocal; c++) {
    if (kind[c] != interior) {
      f(0, ys + c, 0) = 0.0;
    }
  }
}

std::vector<int> Unstructured::wallFacesInBoxes(const std::string& boxes) const {
  std::vector<int> result;
  // Note: ';' starts a comment in BOUT++ input files, so boxes are separated by '|'
  for (const auto& box : strsplit(boxes, '|')) {
    if (trim(box).empty()) {
      continue;
    }
    const auto parts = strsplit(box, ',');
    const std::vector<std::string> numbers(parts.begin(), parts.end());
    if (numbers.size() != 4) {
      throw BoutException(
          "Wall box '{}' should be four numbers in metres: Rmin, Rmax, Zmin, Zmax", box);
    }
    const BoutReal Rmin = stringToReal(trim(numbers[0]));
    const BoutReal Rmax = stringToReal(trim(numbers[1]));
    const BoutReal Zmin = stringToReal(trim(numbers[2]));
    const BoutReal Zmax = stringToReal(trim(numbers[3]));
    for (int i = 0; i < nfaces; i++) {
      if (face_kind[i] == wall and face_R[i] >= Rmin and face_R[i] <= Rmax
          and face_Z[i] >= Zmin and face_Z[i] <= Zmax
          and std::find(result.begin(), result.end(), i) == result.end()) {
        result.push_back(i);
      }
    }
  }
  return result;
}
