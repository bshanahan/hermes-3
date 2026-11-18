#pragma once
#ifndef SHEATH_BOUNDARY_UNSTRUCTURED_H
#define SHEATH_BOUNDARY_UNSTRUCTURED_H

#include "component.hxx"

/// Sheath boundary condition on the wall faces of a flat-cell
/// (unstructured) grid. The physics follows SheathBoundarySimple.
///
/// Every wall face is treated as a sheath, weighted by the magnetic
/// flux through it: a face lying along the field carries no parallel
/// flow. The ghost cell of each wall face is set to the sheath-edge
/// density, temperature, pressure and velocity, and the heat
/// transmitted in excess of the advected energy is removed from the
/// neighbouring cell.
struct SheathBoundaryUnstructured : public NamedComponent<SheathBoundaryUnstructured> {
  SheathBoundaryUnstructured(std::string name, Options& options, Solver*);
  void outputVars(Options& state) override;
  static constexpr auto type = "sheath_boundary_unstructured";

private:
  BoutReal Ge;                    ///< Secondary electron emission coefficient
  BoutReal gamma_e;               ///< Electron sheath heat transmission
  BoutReal gamma_i;               ///< Ion sheath heat transmission
  BoutReal sheath_ion_polytropic; ///< Polytropic coefficient in sheath velocity
  BoutReal density_boundary_mode, pressure_boundary_mode,
      temperature_boundary_mode; ///< BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree
  bool diagnose;                 ///< Save diagnostic variables?

  Field3D hflux_e;              ///< Electron heat loss to the sheath, per volume
  Field3D hflux_i;              ///< Ion heat loss to the sheath (all ions), per volume
  Field3D particle_flux;        ///< Ion particle loss to the sheath (all ions), per volume

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SheathBoundaryUnstructured> registercomponentsheathboundaryunstructured;
}

#endif // SHEATH_BOUNDARY_UNSTRUCTURED_H
