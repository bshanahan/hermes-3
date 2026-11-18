#include "../include/sheath_boundary_unstructured.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/unstructured.hxx"

#include <bout/constants.hxx>
#include <bout/field.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
BoutReal limitFree(BoutReal fm, BoutReal fc, BoutReal mode) {
  if ((fm < fc) && (mode == 0)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  if ((mode == 0) || (mode == 1)) {
    return SQ(fc) / fm; // Exponential
  }
  if (mode == 2) {
    return 2.0 * fc - fm; // Linear
  }
  throw BoutException("Unknown boundary mode");
}

/// Value of `f` at a wall face: midway between the cell and an
/// extrapolation one cell beyond it. `up` is the cell upstream, or -1.
BoutReal sheathValue(const Unstructured::Buffer& f, int up, int cell, BoutReal mode) {
  const BoutReal fc = f[cell];
  if (up < 0) {
    return fc;
  }
  return 0.5 * (limitFree(f[up], fc, mode) + fc);
}
} // namespace

SheathBoundaryUnstructured::SheathBoundaryUnstructured(std::string name,
                                                       Options& alloptions, Solver*)
    : NamedComponent(name, {
                               readIfSet("species:e:AA"),
                               writeBoundaryFinal("species:e:{e_boundary}"),
                               readWrite("species:e:energy_source"),
                               writeBoundaryFinalIfSet("species:e:{e_optional}"),
                               writeBoundaryReadInteriorIfSet("species:e:pressure"),
                               readIfSet("species:{all_species}:charge"),
                               readOnly("species:{ions}:AA"),
                               readWrite("species:{ions}:energy_source"),
                               writeBoundaryFinal("species:{ions}:{ion_boundary}"),
                               writeBoundaryReadInteriorIfSet("species:{ions}:pressure"),
                               writeBoundaryFinalIfSet("species:{ions}:{ion_optional}"),
                           }) {
  Options& options = alloptions[name];

  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);
  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }
  gamma_e = options["gamma_e"]
                .doc("Electron sheath heat transmission coefficient")
                .withDefault(3.5);
  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);
  sheath_ion_polytropic = options["sheath_ion_polytropic"]
                              .doc("Ion polytropic coefficient in Bohm sound speed")
                              .withDefault(1.0);
  density_boundary_mode =
      options["density_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(0);
  pressure_boundary_mode =
      options["pressure_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(0);
  temperature_boundary_mode =
      options["temperature_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(0);
  diagnose = options["diagnose"]
                 .doc("Save additional output diagnostics")
                 .withDefault<bool>(false);

  substitutePermissions("e_boundary", {"density", "temperature"});
  substitutePermissions("e_optional", {"velocity", "momentum"});
  substitutePermissions("ion_boundary", {"density", "temperature"});
  substitutePermissions("ion_optional", {"velocity", "momentum"});
}

void SheathBoundaryUnstructured::transform_impl(GuardedOptions& state) {
  const auto& us = Unstructured::get();
  const int ys = us.ys;

  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];

  Field3D Ne = floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0);
  Field3D Te = GET_NOBOUNDARY(Field3D, electrons["temperature"]);
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
                   ? getNoBoundary<Field3D>(electrons["pressure"])
                   : Te * Ne;
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;
  Field3D Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
                   ? getNoBoundary<Field3D>(electrons["velocity"])
                   : zeroFrom(Ne);
  Field3D NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
                    ? getNoBoundary<Field3D>(electrons["momentum"])
                    : zeroFrom(Ne);
  Ne.allocate();
  Te.allocate();
  Pe.allocate();
  Ve.allocate();
  NVe.allocate();

  // Values on all local cells: the upstream cell may be on another processor
  const auto Ne_b = us.gather(Ne);
  const auto Te_b = us.gather(Te);
  const auto Pe_b = us.gather(Pe);

  // Wall faces
  std::vector<int> wall;
  for (int f = 0; f < us.nfaces; f++) {
    if (us.face_kind[f] == Unstructured::wall) {
      wall.push_back(f);
    }
  }

  // Electron values at the sheath edge of each wall face
  std::vector<BoutReal> nesheath(us.nfaces), tesheath(us.nfaces);
  for (const int f : wall) {
    const int i = us.faceL[f];
    nesheath[f] = sheathValue(Ne_b, us.face_up[f], i, density_boundary_mode);
    tesheath[f] =
        floor(sheathValue(Te_b, us.face_up[f], i, temperature_boundary_mode), 1e-5);
  }

  //////////////////////////////////////////////////////////////////
  // Ion flux into each wall face, to set the electron flux (no net current)

  std::vector<BoutReal> ion_sum(us.nfaces, 0.0);
  for (auto& kv : allspecies.getChildren()) {
    const GuardedOptions species = kv.second;
    if ((kv.first == "e") or !species.isSet("charge")
        or (get<BoutReal>(species["charge"]) == 0.0)) {
      continue; // Skip electrons and non-charged ions
    }
    const Field3D Ni = getNoBoundary<Field3D>(species["density"]);
    const Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
    const BoutReal Mi = getNoBoundary<BoutReal>(species["AA"]);
    const BoutReal Zi = getNoBoundary<BoutReal>(species["charge"]);
    const Field3D Vi = species.isSet("velocity")
                           ? getNoBoundary<Field3D>(species["velocity"])
                           : zeroFrom(Ni);
    const auto Ni_b = us.gather(Ni);
    const auto Ti_b = us.gather(Ti);
    for (const int f : wall) {
      const int i = us.faceL[f];
      const BoutReal sign = (us.flux_over_B[f] >= 0.0) ? 1.0 : -1.0; // b out of the domain?
      const BoutReal nisheath =
          sheathValue(Ni_b, us.face_up[f], i, density_boundary_mode);
      const BoutReal tisheath =
          floor(sheathValue(Ti_b, us.face_up[f], i, temperature_boundary_mode), 1e-5);
      const BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath[f]) / Mi;
      const BoutReal visheath = std::max(sign * Vi(0, ys + i, 0), sqrt(C_i_sq));
      ion_sum[f] += Zi * nisheath * visheath;
    }
  }

  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source = electrons.isSet("energy_source")
                                       ? getNonFinal<Field3D>(electrons["energy_source"])
                                       : zeroFrom(Ne);
  electron_energy_source.allocate();
  hflux_e = zeroFrom(Ne);

  for (const int f : wall) {
    const int i = us.faceL[f];
    const int g = us.faceR[f];
    const BoutReal sign = (us.flux_over_B[f] >= 0.0) ? 1.0 : -1.0;
    const BoutReal area = std::abs(us.flux_over_B[f]); // Area perpendicular to B
    const BoutReal ne = nesheath[f];
    const BoutReal te = tesheath[f];

    // Sheath potential relative to the wall, from zero net current
    const BoutReal phi = te
                         * log(sqrt(te / (Me * TWOPI)) * (1. - Ge) * floor(ne, 1e-5)
                               / floor(ion_sum[f], 1e-5));
    const BoutReal phisheath = floor(phi, 0.0); // Electron saturation at phi = phi_wall
    const BoutReal vesheath = sqrt(te / (TWOPI * Me)) * (1. - Ge) * exp(-phisheath / te);

    // Heat transmitted, less the part carried by advection and kinetic energy
    BoutReal q = gamma_e * te * ne * vesheath;
    hflux_e(0, ys + i, 0) -= q * area / us.vol[i];
    q -= (2.5 * te + 0.5 * Me * SQ(vesheath)) * ne * vesheath;
    electron_energy_source(0, ys + i, 0) -= q * area / us.vol[i];

    Ne(0, ys + g, 0) = ne;
    Te(0, ys + g, 0) = te;
    Pe(0, ys + g, 0) = sheathValue(Pe_b, us.face_up[f], i, pressure_boundary_mode);
    Ve(0, ys + g, 0) = sign * vesheath;
    NVe(0, ys + g, 0) = Me * ne * sign * vesheath;
  }

  setBoundary(electrons["density"], Ne);
  setBoundary(electrons["temperature"], Te);
  setBoundary(electrons["pressure"], Pe);
  set(electrons["energy_source"], electron_energy_source);
  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], Ve);
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], NVe);
  }

  //////////////////////////////////////////////////////////////////
  // Ions

  hflux_i = zeroFrom(Ne);
  particle_flux = zeroFrom(Ne);

  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const
    const BoutReal Zi = species.isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;
    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }
    const BoutReal Mi = get<BoutReal>(species["AA"]);

    Field3D Ni = floor(getNoBoundary<Field3D>(species["density"]), 0.0);
    Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
    Field3D Pi =
        species.isSet("pressure") ? getNoBoundary<Field3D>(species["pressure"]) : Ni * Ti;
    Field3D Vi = species.isSet("velocity") ? getNoBoundary<Field3D>(species["velocity"])
                                           : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum") ? getNoBoundary<Field3D>(species["momentum"])
                                            : Mi * Ni * Vi;
    Field3D energy_source = species.isSet("energy_source")
                                ? getNonFinal<Field3D>(species["energy_source"])
                                : zeroFrom(Ni);
    Ni.allocate();
    Ti.allocate();
    Pi.allocate();
    Vi.allocate();
    NVi.allocate();
    energy_source.allocate();
    const auto Ni_b = us.gather(Ni);
    const auto Ti_b = us.gather(Ti);
    const auto Pi_b = us.gather(Pi);

    for (const int f : wall) {
      const int i = us.faceL[f];
      const int g = us.faceR[f];
      const BoutReal sign = (us.flux_over_B[f] >= 0.0) ? 1.0 : -1.0;
      const BoutReal area = std::abs(us.flux_over_B[f]);

      const BoutReal nisheath =
          sheathValue(Ni_b, us.face_up[f], i, density_boundary_mode);
      const BoutReal tisheath =
          floor(sheathValue(Ti_b, us.face_up[f], i, temperature_boundary_mode), 1e-5);
      const BoutReal pisheath =
          sheathValue(Pi_b, us.face_up[f], i, pressure_boundary_mode);

      // Ion speed into the sheath: at least the sound speed
      const BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath[f]) / Mi;
      const BoutReal visheath = std::max(sign * Vi(0, ys + i, 0), sqrt(C_i_sq));

      BoutReal q = gamma_i * tisheath * nisheath * visheath;
      hflux_i(0, ys + i, 0) -= q * area / us.vol[i];
      particle_flux(0, ys + i, 0) -= nisheath * visheath * area / us.vol[i];
      q -= (2.5 * tisheath + 0.5 * Mi * SQ(visheath)) * nisheath * visheath;
      energy_source(0, ys + i, 0) -= q * area / us.vol[i];

      Ni(0, ys + g, 0) = nisheath;
      Ti(0, ys + g, 0) = tisheath;
      Pi(0, ys + g, 0) = pisheath;
      Vi(0, ys + g, 0) = sign * visheath;
      NVi(0, ys + g, 0) = Mi * nisheath * sign * visheath;
    }

    setBoundary(species["density"], Ni);
    setBoundary(species["temperature"], Ti);
    setBoundary(species["pressure"], Pi);
    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], Vi);
    }
    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], NVi);
    }
    set(species["energy_source"], energy_source);
  }
}

void SheathBoundaryUnstructured::outputVars(Options& state) {
  if (!diagnose) {
    return;
  }
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm;

  set_with_attrs(state["Ee_sheath"], hflux_e,
                 {{"time_dimension", "t"},
                  {"units", "W m^-3"},
                  {"conversion", Pnorm * Omega_ci},
                  {"long_name", "electron energy loss to the sheath"},
                  {"source", "sheath_boundary_unstructured"}});
  set_with_attrs(state["Ei_sheath"], hflux_i,
                 {{"time_dimension", "t"},
                  {"units", "W m^-3"},
                  {"conversion", Pnorm * Omega_ci},
                  {"long_name", "ion energy loss to the sheath, all ions"},
                  {"source", "sheath_boundary_unstructured"}});
  set_with_attrs(state["Si_sheath"], particle_flux,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"long_name", "ion particle loss to the sheath, all ions"},
                  {"source", "sheath_boundary_unstructured"}});
}
