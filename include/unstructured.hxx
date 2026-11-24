#pragma once
#ifndef HERMES_UNSTRUCTURED_HXX
#define HERMES_UNSTRUCTURED_HXX

#include <bout/field3d.hxx>
#include <bout/mesh.hxx>
#include <bout/utils.hxx>

#include <algorithm>
#include <cmath>
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
  std::vector<BoutReal> cellR, cellZ; ///< Cell centre (ghost: face mid-point), normalised
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
  std::vector<BoutReal> kpar;        ///< Parallel transmissibility (two-point gradient)
  /// Parallel conduction from the diamond gradient: centre and vertex coefficients
  std::vector<BoutReal> kparn, vpar;
  bool parallel_correction{false};   ///< Use kparn, vpar instead of kpar
  std::vector<BoutReal> tperp;       ///< Radial transmissibility
  std::vector<BoutReal> tiso;        ///< Transmissibility for the whole of grad_perp
  std::vector<BoutReal> face_area;   ///< Surface area per radian
  std::vector<BoutReal> face_length; ///< Poloidal length
  std::vector<BoutReal> nR, nZ;      ///< Unit normal in the poloidal plane, from L to R
  std::vector<BoutReal> face_R, face_Z; ///< Position of the face mid-point [m]
  /// Area of a wall face projected on the flux surface normal, per radian.
  /// Zero on other faces.
  std::vector<BoutReal> wall_perp_area;

  /// Drifts (see flatgrid.py). Empty if the grid file has no drift geometry.
  /// E x B volume flow from L to R = face_drift * (phi at vertex 2 - vertex 1);
  /// flux of curl(b/B) from L to R = face_curl. Both include `field_sign`.
  std::vector<BoutReal> face_drift, face_curl;
  bool has_drifts{false};   ///< Drift geometry and vertex stencils loaded?
  BoutReal field_sign{1.0}; ///< -1 reverses the toroidal field
  bool exb_through_wall{true}; ///< Does the E x B drift carry plasma into the wall?
  /// Drift fluxes with limited second-order face values (otherwise first-order upwind)?
  bool drift_second_order{true};

  /// E x B volume flow through every face, from L to R, for the potential
  /// `phi`. Through wall faces only the flow into the wall is kept
  /// (none at all if `exb_through_wall` is false).
  std::vector<BoutReal> flowExB(const Field3D& phi) const;
  /// Flow of the magnetic drift through every face for a species with
  /// temperature `T` and charge `charge`: T / charge * curl(b/B) . dS.
  /// Zero through wall faces.
  std::vector<BoutReal> flowMagneticDrift(const Field3D& T, BoutReal charge) const;
  /// Area-weighted mean over the core boundary of the ghost-cell values of a
  /// gathered field (all processors)
  BoutReal coreMean(const Buffer& f) const;
  /// Div(f v) for a given volume flow through each face. On interior faces
  /// `f` is taken from the upwind side: with `drift_second_order` a limited
  /// linear reconstruction from the upwind cell (TVD, van Albada limiter),
  /// otherwise the upwind cell value. Core faces carry the mean of the core
  /// boundary values, in either direction; wall faces the value of their cell.
  /// If `boundary_flow` is given, what leaves through each boundary face
  /// (wall or core) is added to the ghost cell of that face.
  /// If `weight` is given (one value per face, 0 to 1), the second-order
  /// face value is blended with the upwind one: weight 0 = upwind.
  Field3D Div_f_flow(const Field3D& f, const std::vector<BoutReal>& flow,
                     Field3D* boundary_flow = nullptr,
                     const std::vector<BoutReal>* weight = nullptr) const;
  /// Weight of the second-order drift face values: n_f / (n_f + n_blend),
  /// with n_f the harmonic mean of the densities `n` of the two cells and
  /// n_blend = drift_second_order_density. In tenuous plasma the drift
  /// fluxes then fall back to upwind values, which keep an implicit step
  /// positive (second-order values do not). With drift_second_order_closed_only
  /// the weight is zero on every face that is not between two cells of
  /// closed flux tubes.
  std::vector<BoutReal> secondOrderWeight(const Field3D& n) const;
  BoutReal drift_second_order_density{0.0}; ///< Normalised; 0 = second order everywhere
  /// Second-order drift face values only between cells of closed flux tubes?
  bool drift_second_order_closed_only{true};
  std::vector<char> closed; ///< Is the cell in a closed flux tube? (local cells, halo included)
  /// Green-Gauss gradient (R, Z components, per normalised length) of a
  /// gathered field in every local cell, halo included (communicates).
  /// Zero in ghost and padding cells. `like` only sets the field's mesh.
  void cellGradient(const Field3D& like, const Buffer& f, Buffer& gR, Buffer& gZ) const;

  /// Skew correction of the cross-field operators (see flatgrid.py): the
  /// gradient at a face also takes the difference of the field between the
  /// two vertices of the face, interpolated from nearby cells.
  bool skew_correction{false};
  bool skew_correction_mean{false}; ///< Also in Div_a_Grad_perp, which the neutrals use
  std::vector<BoutReal> vperp, viso; ///< Coefficient of the vertex difference
  std::vector<int> sptr, scell;      ///< Vertex difference as a sum over cells (CSR by face)
  std::vector<BoutReal> sw;
  /// Two-point gradient term plus the bounded vertex correction.
  /// The correction C is kept below
  ///     M = skew_limit * sqrt(A^2 + (skew_smooth * trans * h)^2),
  /// A the two-point term and h a mean of |f| in the two cells. The term in h
  /// lets a resolved profile (cells smaller than about skew_smooth times its
  /// gradient length) keep its whole correction, so that the operator stays
  /// second order; the bound acts only where the profile is under-resolved.
  /// With the harmonic mean h -> 0 when either cell empties, and then |C| < |A|:
  /// the flow into an empty cell never reverses (positivity). The bound is
  /// approached smoothly, C / (1 + |C/M|^p)^(1/p), which Newton iterations need.
  /// skew_smooth = 0.05, arithmetic mean, p = 2 is the bound of the first version.
  /// The size h of the field that sets the floor of the bound depends on what it is:
  enum class Scale {
    positive,  ///< density, temperature, pressure: harmonic mean of the two cells
    logarithm, ///< a logarithm: h = 1, a change by a factor e (normalisation independent)
    signed_    ///< velocity, potential: root-mean-square of the two cells (smooth through 0)
  };
  BoutReal skew_limit{1.0};
  BoutReal skew_smooth{1.0};
  bool skew_harmonic{true}; ///< false: arithmetic mean of |f| for every field (first version)
  BoutReal skew_power{2.0};
  BoutReal correctedGradient(BoutReal two_point, BoutReal trans, BoutReal v,
                             const Buffer& f, int i, Scale scale = Scale::positive) const {
    const BoutReal fa = f[faceL[i]];
    const BoutReal fb = f[faceR[i]];
    const BoutReal a = std::abs(fa);
    const BoutReal b = std::abs(fb);
    BoutReal h = 0.5 * (a + b);
    if (skew_harmonic) {
      h = (scale == Scale::positive)    ? 2.0 * a * b / (a + b + 1e-300)
          : (scale == Scale::logarithm) ? 1.0
                                        : std::sqrt(0.5 * (SQ(fa) + SQ(fb)));
    }
    const BoutReal bound =
        skew_limit * std::sqrt(SQ(two_point) + SQ(skew_smooth * h * trans));
    const BoutReal correction = v * vertexDifference(f, i);
    if (skew_power == 2.0) {
      return two_point
             + correction * bound / std::sqrt(SQ(bound) + SQ(correction) + 1e-300);
    }
    const BoutReal r = std::abs(correction) / (bound + 1e-300);
    return two_point + correction / std::pow(1.0 + std::pow(r, skew_power), 1.0 / skew_power);
  }
  /// Value at vertex 2 minus value at vertex 1 of face `i`
  BoutReal vertexDifference(const Buffer& f, int i) const {
    BoutReal sum = 0.0;
    for (int k = sptr[i]; k < sptr[i + 1]; k++) {
      sum += sw[k] * f[scell[k]];
    }
    return sum;
  }

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
  /// With the skew correction the sign of the whole gradient picks the side.
  /// If `boundary_flow` is given, the flow out of the plasma through each
  /// core face is added to its ghost cell. `scale` says what kind of field f is,
  /// for the bound on the skew correction.
  Field3D Div_a_Grad_perp_upwind(const Field3D& a, const Field3D& f,
                                 bool isotropic = false,
                                 Field3D* boundary_flow = nullptr,
                                 Scale scale = Scale::positive) const;
  /// Div(a Grad_perp(f)) with the mean of `a` at each face. Includes the
  /// poloidal part of Grad_perp. No flux through wall faces.
  Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& f,
                          Scale scale = Scale::positive) const;
  /// Loss through wall faces with a decay length: the flow out of each cell
  /// next to the wall is a * f * (wall area projected on the flux surface
  /// normal) * inv_length, with `a` and `f` taken in that cell. Returns the
  /// (negative) source. If `flow` is given, the flow through each wall face
  /// is added to its ghost cell.
  Field3D wallLoss(const Field3D& a, const Field3D& f, BoutReal inv_length,
                   Field3D* flow = nullptr) const;
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
