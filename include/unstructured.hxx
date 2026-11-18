#pragma once
#ifndef HERMES_UNSTRUCTURED_HXX
#define HERMES_UNSTRUCTURED_HXX

#include <bout/field3d.hxx>
#include <bout/mesh.hxx>

#include <optional>
#include <string>
#include <utility>
#include <vector>

/// Flat-cell ("unstructured") finite-volume grid and its operators.
///
/// The cells are stored along the BOUT++ y index of a mesh with nx = 1,
/// nz = 1 and no x guard cells. Each real cell is followed by one ghost
/// cell per boundary face it has; a ghost cell holds the value AT its
/// boundary face. Inert padding cells make the total a multiple of the
/// number of processors.
///
/// The grid file is written by tools/unstructured/flatgrid.py, which
/// documents the geometric quantities. Everything is normalised on loading.
///
/// In parallel each processor owns a contiguous piece of the list. Cells
/// are numbered locally: owned cells first (cell `c` is at y index
/// `ys + c`), then the halo cells owned by other processors that the
/// operators need. A `Buffer` holds one value per local cell, halo
/// included; `gather` fills it from a field. A cell and its ghost cells
/// are always owned by the same processor.
///
/// Operators only write to real cells; other cells of the result are zero.
class Unstructured {
public:
  enum Kind { interior = 0, core = 1, wall = 2, padding = 3 };
  using Buffer = std::vector<BoutReal>;

  /// Does the mesh file hold a flat-cell grid?
  static bool present(Mesh* mesh);
  /// Load the grid, and set the cell volumes and field strength in the
  /// BOUT++ coordinates. Lengths are divided by rho_s0, fields by Bnorm.
  static void initialise(Mesh* mesh, BoutReal rho_s0, BoutReal Bnorm);
  /// Has initialise() been called?
  static bool active() { return instance.has_value(); }
  static const Unstructured& get();

  int nlocal{0}; ///< Number of cells owned by this processor
  int nall{0};   ///< Owned cells plus halo cells
  int nfaces{0}; ///< Faces with at least one owned cell
  int ys{0};     ///< y index of owned cell 0
  int global_start{0}; ///< Position of owned cell 0 in the whole list
  bool second_order{true}; ///< Limited linear reconstruction along flux tubes?
  BoutReal length_norm{1.0}; ///< rho_s0 [m]

  // Cell data, length nall
  std::vector<BoutReal> vol;    ///< Volume per radian
  std::vector<BoutReal> B;      ///< Field strength
  std::vector<BoutReal> area;   ///< Poloidal area
  std::vector<BoutReal> bR, bZ; ///< Poloidal components of the unit vector along B
  std::vector<int> prev, next;  ///< Neighbours along the flux tube, or -1
  /// Ghost cell closing the flux tube beyond a cell that has no `prev` / `next`
  /// neighbour, or -1. Used like a guard cell when reconstructing face values.
  std::vector<int> ghost_prev, ghost_next;
  /// interior for a real cell; core or wall for a ghost cell on that boundary; padding
  std::vector<int> kind;

  // Face data, length nfaces. Boundary faces always have both cells owned.
  std::vector<int> faceL, faceR;     ///< Cells either side. R is a ghost on boundary faces
  std::vector<int> face_kind;        ///< interior / core / wall
  std::vector<int> sideL, sideR;     ///< +1 if the other cell is `next`, -1 if `prev`, else 0
  std::vector<int> face_up;          ///< Boundary faces: cell upstream of L, or -1
  std::vector<BoutReal> flux_over_B; ///< Magnetic flux from L to R, divided by |B|
  std::vector<BoutReal> kpar;        ///< Parallel transmissibility
  std::vector<BoutReal> tperp;       ///< Radial transmissibility
  std::vector<BoutReal> tiso;        ///< Transmissibility for the whole of grad_perp
  std::vector<BoutReal> face_area;   ///< Surface area per radian
  std::vector<BoutReal> face_length; ///< Poloidal length
  std::vector<BoutReal> nR, nZ;      ///< Unit normal in the poloidal plane, from L to R
  std::vector<BoutReal> face_R, face_Z; ///< Position of the face mid-point [m]

  /// Values of `f` on all local cells, halo cells included (communicates)
  Buffer gather(const Field3D& f) const;

  /// Divergence of a parallel flow, Div(b f v), with a Lax flux.
  /// Boundary faces carry the flux given by the ghost values.
  Field3D Div_par_mod(const Field3D& f, const Field3D& v, const Field3D& wave_speed) const;
  /// Div(b f v v), for the parallel momentum equation.
  /// With `fixflux` the flux through boundary faces is that of the ghost
  /// values; otherwise it is built from the cell's own face values, with a
  /// penalty on the velocity difference (as FV::Div_par_fvv does).
  Field3D Div_par_fvv(const Field3D& f, const Field3D& v, const Field3D& wave_speed,
                      bool fixflux = true) const;
  /// Div(b f), central
  Field3D Div_par(const Field3D& f) const;
  /// b . Grad(f), central
  Field3D Grad_par(const Field3D& f) const;
  /// Div(b K b.Grad(f)). No flux through boundary faces.
  Field3D Div_par_K_Grad_par(const Field3D& K, const Field3D& f) const;
  /// Div(a Grad_perp(f)), with `a` taken from the upwind side.
  /// Radial only by default; `isotropic` includes the poloidal part of Grad_perp.
  /// No flux through wall faces; core faces use the ghost value.
  Field3D Div_a_Grad_perp_upwind(const Field3D& a, const Field3D& f,
                                 bool isotropic = false) const;
  /// Div(a Grad_perp(f)) with the mean of `a` at each face. Includes the
  /// poloidal part of Grad_perp. No flux through wall faces.
  Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& f) const;
  /// |Grad_perp(f)|^2 at cell centres, from a Green-Gauss gradient
  Field3D gradPerpSq(const Field3D& f) const;

  /// Values of a gathered field either side of owned cell `c` along its
  /// flux tube, in the role of f[i.ym()] and f[i.yp()] on a structured grid.
  /// At the end of a tube the ghost cell value is reflected, like a guard
  /// cell; without a neighbour the cell's own value is returned.
  std::pair<BoutReal, BoutReal> tubeValues(const Buffer& buf, int c) const {
    const BoutReal fc = buf[c];
    if (kind[c] != interior) {
      return {fc, fc};
    }
    const BoutReal fm = (prev[c] >= 0)
                            ? buf[prev[c]]
                            : ((ghost_prev[c] >= 0) ? 2. * buf[ghost_prev[c]] - fc : fc);
    const BoutReal fp = (next[c] >= 0)
                            ? buf[next[c]]
                            : ((ghost_next[c] >= 0) ? 2. * buf[ghost_next[c]] - fc : fc);
    return {fm, fp};
  }

  /// Set ghost cells: copy the neighbouring cell (zero gradient), except
  /// on the core boundary if a value is given.
  void setGhosts(Field3D& f, std::optional<BoutReal> core_value = std::nullopt) const;
  /// Set every cell that is not a real cell (ghost and padding cells) to zero
  void zeroGhosts(Field3D& f) const;

  /// Wall faces (local indices) whose mid-point lies in one of the boxes
  /// "Rmin,Rmax,Zmin,Zmax | ..." given in metres
  std::vector<int> wallFacesInBoxes(const std::string& boxes) const;

private:
  static std::optional<Unstructured> instance;

  /// Halo exchange: for each other processor, the owned cells it needs and
  /// where its cells go in a Buffer
  struct Exchange {
    int rank;
    std::vector<int> send; ///< Owned cells (local index) to send
    int recv_start;        ///< First local index of the cells received
    int recv_count;
  };
  std::vector<Exchange> exchanges;

  /// Values of a gathered field at the two tube-direction faces of each cell
  void reconstruct(const Buffer& f, Buffer& at_prev, Buffer& at_next) const;
};

#endif // HERMES_UNSTRUCTURED_HXX
