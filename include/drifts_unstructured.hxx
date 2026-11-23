#pragma once
#ifndef DRIFTS_UNSTRUCTURED_H
#define DRIFTS_UNSTRUCTURED_H

#include "component.hxx"

#include <vector>

/// Drifts and electrostatic potential on a flat-cell (unstructured) grid,
/// for axisymmetric transport calculations.
///
/// The potential phi (zero at the wall) is evolved in pseudo-time until the
/// current is divergence free:
///
///   d phi / dt = -lambda * Div(J),
///   J = b J_par + J_drift + J_perp
///
///   J_par   = sigma_par * (-d_par phi + d_par p_e / n_e + 0.71 d_par T_e)
///             Ohm's law for massless electrons, sigma_par = n_e / (0.51 m_e nu_ei).
///             Through wall faces: the sheath current of all charged species.
///   J_drift = sum over species of p_s * curl(b / B): the current of the
///             magnetic (grad-B and curvature) drifts
///   J_perp  = -sigma_perp * (grad_perp phi + grad_perp p_i / n_e),
///             sigma_perp = e * n_e * perp_mobility. An anomalous current that
///             damps the perpendicular ion flow (E x B plus diamagnetic). It
///             makes the potential on closed flux surfaces well defined.
///
/// Every charged species is advected by the E x B drift b x grad(phi) / B and
/// by its magnetic drift (T / q) curl(b / B), both through the faces of the
/// poloidal cells (see tools/unstructured/flatgrid.py):
///
///   density     -Div(n (v_E + v_D))
///   pressure    (3/2) d p / dt = -Div((3/2) p v_E) - Div((5/2) p v_D)
///   momentum    -Div(m n v_par (v_E + v_D))
///
/// as in the structured-grid components (evolve_*, diamagnetic_drift).
/// The E x B drift carries plasma into the wall where it points that way
/// (unstructured:exb_through_wall), in addition to the flow along the field
/// through the sheath; those ions are recycled and sputter like the others
/// (species:<ion>:wall_drift_flow). Nothing drifts out of the wall, and the
/// magnetic drift does not cross it. The polarisation drift is not included.
///
/// The electron parallel velocity is set to that of the ion flow, so this
/// component replaces zero_current. What the current changes in the electron
/// energy balance is added from the current at the cell faces:
/// Div((5/2 + 0.71) T_e J_par b) + J_par E_par. It must come after the species and
/// before the sheath boundary in the list of components. The parallel force
/// on the ions is still taken from electron_force_balance.
///
/// Options, in the section named after the component:
///   lambda          Relaxation rate of the potential (normalised)
///   perp_mobility   Anomalous perpendicular mobility [m^2 / V / s]
///   ion_pressure_term     Include grad p_i / n_e in J_perp? (default true)
///   exb_advection         E x B advection? (default true)
///   magnetic_drift        Magnetic drift advection and current? (default true)
///   current_in_electron_energy     Heat carried by the current, and J.E heating? (default true)
///   diagnose
/// The sign of the toroidal field is unstructured:field_sign.
struct DriftsUnstructured : public NamedComponent<DriftsUnstructured> {
  DriftsUnstructured(std::string name, Options& alloptions, Solver* solver);
  void finally(const Options& state) override;
  void outputVars(Options& state) override;
  static constexpr auto type = "drifts_unstructured";

private:
  Field3D phi; ///< Potential, evolving. Normalised to Tnorm

  BoutReal lambda;        ///< Relaxation rate
  BoutReal perp_mobility; ///< Normalised: sigma_perp = n_e * perp_mobility
  BoutReal conductivity_factor; ///< Multiplies the parallel conductivity (for tests)
  bool ion_pressure_term, exb_advection, magnetic_drift, electron_current, diagnose;
  BoutReal exb_wall_density; ///< Normalised: E x B wall outflow fades below this density
  bool second_order_exb; ///< Second-order face values for E x B as well as the magnetic drift?
  BoutReal drift_density; ///< Normalised: drift flows fade out below this density
  BoutReal Tnorm, Nnorm, Omega_ci;

  Field3D sigma_par;  ///< Parallel conductivity
  Field3D Jpar;       ///< Parallel current at cell centres
  Field3D Ve, NVe;    ///< Electron parallel velocity and momentum
  Field3D DivJdrift;  ///< Divergence of the magnetic drift current
  Field3D Pi_total;   ///< Pressure of all positive ions
  Field3D DivJpar, DivJperp, DivJsheath; ///< Diagnostics
  Field3D Sn_exb, Sn_mag; ///< Diagnostics: electron particle sources of the drifts
  Field3D wall_particles; ///< Ions carried into the wall by the E x B drift, per volume
  Field3D wall_energy;    ///< Thermal energy carried into the wall by it, per volume

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<DriftsUnstructured> registercomponentdriftsunstructured;
}

#endif // DRIFTS_UNSTRUCTURED_H
