#pragma once
#ifndef BACKGROUND_PHI_H
#define BACKGROUND_PHI_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Adds terms to Density, Pressure and Momentum equations following the
/// stellarator 2-point model from Yuhe Feng et al., PPCF 53 (2011) 024009
/// The terms add the effective parallel contribution of perpendicular drifts
/// in long connection length scenarios.
/// B Shanahan 2023 <brendan.shanahan@ipp.mpg.de>
struct BackgroundPhi : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - D           Perpendicular density diffusion coefficient
  ///   - chi         Perpendicular heat diffusion coefficient
  ///   - nu          Perpendicular momentum diffusion coefficient
  ///   - Vbn         Binormal velocity from cross-field drifts
  ///   - Theta       Field line pitch as described by Feng et al.
  ///
  BackgroundPhi(std::string name, Options& options, Solver* solver);
  
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure correction 
  ///     - momentum correction
  ///     - density correction
  ///
  void transform(Options& state) override;


private:
  std::string name; ///< Short name of the species e.g. h+
  Field3D phi_bg; ///< background potential from imposed electric field.

};

namespace {
RegisterComponent<BackgroundPhi> registercomponentbackgroundphi("background_phi");
}

#endif // BACKGROUND_PHI
