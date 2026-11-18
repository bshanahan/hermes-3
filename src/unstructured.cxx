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

  std::vector<std::vector<int>> adjacent(ncells);
  for (int f = 0; f < nfaces; f++) {
    adjacent[all_faceL[f]].push_back(all_faceR[f]);
    adjacent[all_faceR[f]].push_back(all_faceL[f]);
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
  us.kind = pick(all_kind, cells);
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
  us.faceL.resize(us.nfaces);
  us.faceR.resize(us.nfaces);
  us.face_up.resize(us.nfaces);
  for (int i = 0; i < us.nfaces; i++) {
    const int f = faces[i];
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
               "{} cells, {} halo cells and {} faces on this processor\n",
               nreal, nghost, nfaces, us.nlocal, halo.size(), us.nfaces);
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
                                             bool isotropic) const {
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
    const BoutReal gradient = trans[i] * (f[R] - f[L]);
    // Flow is from high to low f: coefficient from the cell it leaves
    const BoutReal flow = gradient * ((gradient > 0) ? a[R] : a[L]); // Into L
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (face_kind[i] == interior and R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
    }
  }
  return result;
}

Field3D Unstructured::Div_a_Grad_perp(const Field3D& a_in, const Field3D& f_in) const {
  Field3D result{zeroFrom(f_in)};
  const Buffer a = gather(a_in);
  const Buffer f = gather(f_in);
  for (int i = 0; i < nfaces; i++) {
    if (tiso[i] == 0.0) {
      continue;
    }
    const int L = faceL[i];
    const int R = faceR[i];
    const BoutReal flow = 0.5 * (a[L] + a[R]) * tiso[i] * (f[R] - f[L]); // Into L
    if (L < nlocal) {
      result(0, ys + L, 0) += flow / vol[L];
    }
    if (face_kind[i] == interior and R < nlocal) {
      result(0, ys + R, 0) -= flow / vol[R];
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
