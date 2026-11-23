#include "../include/drifts_unstructured.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/unstructured.hxx"

#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/solver.hxx>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

DriftsUnstructured::DriftsUnstructured(std::string name, Options& alloptions,
                                       Solver* solver)
    : NamedComponent(name, {readIfSet("species:{all_species}:charge"),
                            readIfSet("species:{all_species}:AA"),
                            readIfSet("species:{all_species}:{inputs}", Regions::Interior),
                            readWrite("species:{all_species}:{outputs}"),
                            readWrite("species:e:velocity"),
                            readWrite("species:e:momentum"),
                            readWrite("fields:phi_unstructured")}) {
  if (!Unstructured::active()) {
    throw BoutException("drifts_unstructured needs a flat-cell (unstructured) grid");
  }
  if (!Unstructured::get().has_drifts) {
    throw BoutException("drifts_unstructured: the grid file has no drift geometry "
                        "(us_face_drift, us_face_curl, vertex stencils). Regenerate it");
  }
  Options& options = alloptions[name];
  const auto& units = alloptions["units"];
  Tnorm = get<BoutReal>(units["eV"]);
  Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);
  const BoutReal Bnorm = get<BoutReal>(units["Tesla"]);

  lambda = options["lambda"]
               .doc("Relaxation rate of the potential: d phi/dt = -lambda * Div(J), "
                    "in normalised units. The steady state does not depend on it")
               .withDefault<BoutReal>(1e3);
  // sigma_perp = e n mu, normalised to e Nnorm / Bnorm
  perp_mobility = options["perp_mobility"]
                      .doc("Anomalous perpendicular mobility [m^2/V/s]: "
                           "J_perp = -e n_e mobility (grad_perp phi + grad_perp p_i / n_e)")
                      .withDefault<BoutReal>(1e-4)
                  * Bnorm;
  conductivity_factor =
      options["conductivity_factor"]
          .doc("Multiplies the parallel conductivity. For tests: with a very small "
               "value and no perpendicular mobility every sheath floats, and the "
               "solution without drifts must be recovered")
          .withDefault<BoutReal>(1.0);
  ion_pressure_term = options["ion_pressure_term"]
                          .doc("Include the ion pressure gradient in the anomalous "
                               "perpendicular current?")
                          .withDefault<bool>(true);
  drift_density = options["drift_density"]
                      .doc("Density [m^-3] below which the drift flows fade out (factor "
                           "n/(n + this) on every face), for near-vacuum far-SOL cells. "
                           "0: no fading")
                      .withDefault<BoutReal>(1e17)
                  / Nnorm;
  second_order_exb = options["second_order_exb"]
                         .doc("Second-order face values for the E x B fluxes too (with "
                              "unstructured:drift_second_order)? If not, E x B stays "
                              "first-order upwind and only the magnetic drift is second order")
                         .withDefault<bool>(false);
  exb_wall_density = options["exb_wall_density"]
                         .doc("Density [m^-3] below which E x B flow into the wall fades "
                              "out (factor n/(n + this)). 0: no fading")
                         .withDefault<BoutReal>(1e16)
                     / Nnorm;
  exb_advection =
      options["exb_advection"].doc("Include E x B advection?").withDefault<bool>(true);
  magnetic_drift = options["magnetic_drift"]
                       .doc("Include the magnetic (grad-B and curvature) drifts and "
                            "their current?")
                       .withDefault<bool>(true);
  electron_current = options["current_in_electron_energy"]
                         .doc("Include the heat carried by the parallel current and "
                              "the work J_par E_par in the electron energy equation?")
                         .withDefault<bool>(true);
  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  solver->add(phi, "phi");

  substitutePermissions("inputs", {"density", "temperature", "pressure", "velocity"});
  substitutePermissions("outputs", {"density_source", "energy_source", "momentum_source",
                                    "wall_drift_flow", "core_energy_flow"});
}

void DriftsUnstructured::transform_impl(GuardedOptions& state) {
  const auto& us = Unstructured::get();
  us.setGhosts(phi); // Zero gradient: no current through the boundary by conduction

  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];
  const Field3D Ne = getNoBoundary<Field3D>(electrons["density"]);
  const Field3D Te = getNoBoundary<Field3D>(electrons["temperature"]);
  const Field3D Pe = getNoBoundary<Field3D>(electrons["pressure"]);
  const BoutReal Me = get<BoutReal>(electrons["AA"]);

  //////////////////////////////////////////////////////////////////
  // Ion flow, and the friction that sets the parallel conductivity

  Field3D ion_current{zeroFrom(Ne)};  // sum of Z n v
  Field3D friction{zeroFrom(Ne)};     // sum of coefficient * nu_ei, normalised
  Field3D thermal{zeroFrom(Ne)};      // sum of 0.71 Z^2 n_i
  Pi_total = zeroFrom(Ne);
  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first];
    if (kv.first == "e" or !species.isSet("charge")
        or (get<BoutReal>(species["charge"]) <= 0.0)) {
      continue;
    }
    const BoutReal Zi = get<BoutReal>(species["charge"]);
    const BoutReal Ai = get<BoutReal>(species["AA"]);
    const Field3D Ni = getNoBoundary<Field3D>(species["density"]);
    const Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
    if (species.isSet("velocity")) {
      ion_current += Zi * Ni * getNoBoundary<Field3D>(species["velocity"]);
    }
    Pi_total += species.isSet("pressure") ? getNoBoundary<Field3D>(species["pressure"])
                                          : Field3D{Ni * Ti};
    thermal += 0.71 * SQ(Zi) * Ni;
    // Friction coefficient of Braginskii (as braginskii_friction)
    const BoutReal coef = Zi == 1 ? 0.51 : Zi == 2 ? 0.44 : Zi == 3 ? 0.40 : 0.38;
    const BoutReal me_mi = SI::Me / (SI::Mp * Ai);
    for (int c = 0; c < us.nlocal; c++) {
      // Electron-ion collision frequency, as braginskii_collisions
      const BoutReal te = Te(0, us.ys + c, 0) * Tnorm;
      const BoutReal ti = Ti(0, us.ys + c, 0) * Tnorm;
      const BoutReal ne = Ne(0, us.ys + c, 0) * Nnorm;
      const BoutReal ni = Ni(0, us.ys + c, 0) * Nnorm;
      const BoutReal Te_lim = std::max(te, 0.1);
      const BoutReal Ti_lim = std::max(ti, 0.1);
      const BoutReal Ne_lim = std::max(ne, 1e10);
      const BoutReal Ni_lim = std::max(ni, 1e10);
      const BoutReal coulomb_log =
          (Te_lim < Ti_lim * me_mi)
              ? 23 - 0.5 * log(Ni_lim) + 1.5 * log(Ti_lim) - log(SQ(Zi) * Ai)
          : (Te_lim < exp(2) * SQ(Zi))
              ? 30.0 - 0.5 * log(Ne_lim) - log(Zi) + 1.5 * log(Te_lim)
              : 31.0 - 0.5 * log(Ne_lim) + log(Te_lim);
      const BoutReal vesq = 2 * softFloor(te, 0.1) * SI::qe / SI::Me;
      const BoutReal visq = 2 * softFloor(ti, 0.1) * SI::qe / (SI::Mp * Ai);
      const BoutReal nu = SQ(SQ(SI::qe) * Zi) * std::max(ni, 0.0)
                          * softFloor(coulomb_log, 1.0) * (1. + me_mi)
                          / (3 * pow(PI * (vesq + visq), 1.5) * SQ(SI::e0 * SI::Me));
      friction(0, us.ys + c, 0) += coef * nu / Omega_ci;
    }
  }

  //////////////////////////////////////////////////////////////////
  // Ohm's law at the cell centres

  sigma_par = zeroFrom(Ne);
  Field3D thermal_coef{zeroFrom(Ne)};
  for (int c = 0; c < us.nlocal; c++) {
    const int y = us.ys + c;
    const BoutReal ne = softFloor(Ne(0, y, 0), 1e-7);
    sigma_par(0, y, 0) =
        conductivity_factor * ne / (Me * std::max(friction(0, y, 0), 1e-20));
    thermal_coef(0, y, 0) = thermal(0, y, 0) / ne;
  }
  us.setGhosts(sigma_par);
  us.setGhosts(thermal_coef);

  // The electrons move with the ions here. What the current adds to the
  // electron energy balance is put in below, from the current at the faces.
  Ve = zeroFrom(Ne);
  for (int c = 0; c < us.nlocal; c++) {
    const int y = us.ys + c;
    Ve(0, y, 0) = ion_current(0, y, 0) / softFloor(Ne(0, y, 0), 1e-7);
  }
  us.setGhosts(Ve);
  NVe = Me * Ne * Ve;
  set(electrons["velocity"], Ve);
  set(electrons["momentum"], NVe);

  // Current between neighbouring cells along the field,
  //   I = sigma * kpar * (-(phi_R - phi_L) + (p_R - p_L) / n + c_T (T_R - T_L)),
  // its divergence, and the electron heating
  //   Div((5/2 + c_T) T_e J_par b) + J_par E_par
  // (heat carried by the electrons against the current, plus the work done by
  // the field). c_T = 0.71 for singly charged ions.
  {
    const auto phi_b = us.gather(phi);
    const auto pe_b = us.gather(Pe);
    const auto te_b = us.gather(Te);
    const auto ne_b = us.gather(Ne);
    const auto sigma_b = us.gather(sigma_par);
    const auto ct_b = us.gather(thermal_coef);
    DivJpar = zeroFrom(Ne);
    Jpar = zeroFrom(Ne);
    Field3D heating{zeroFrom(Ne)};
    Field3D weight{zeroFrom(Ne)};
    for (int f = 0; f < us.nfaces; f++) {
      if (us.face_kind[f] != Unstructured::interior or us.kpar[f] == 0.0) {
        continue;
      }
      const int L = us.faceL[f];
      const int R = us.faceR[f];
      const BoutReal sigma = 0.5 * (sigma_b[L] + sigma_b[R]);
      const BoutReal ne = std::max(0.5 * (ne_b[L] + ne_b[R]), 1e-7);
      const BoutReal ct = 0.5 * (ct_b[L] + ct_b[R]);
      // Current from L to R
      const BoutReal current =
          sigma * us.kpar[f]
          * ((phi_b[L] - phi_b[R]) - (pe_b[L] - pe_b[R]) / ne - ct * (te_b[L] - te_b[R]));
      // Electrons run against the current and bring their heat with them
      // The energy the current gives or takes from the electrons fades out in
      // near-vacuum plasma, as the drifts do: there the work J.E does not scale
      // with the (tiny) electron energy and drives P_e negative within a
      // Newton iteration. The current itself is kept, so charge is conserved.
      BoutReal fade_e = 1.0;
      if (drift_density > 0.0) {
        const BoutReal nmin = 2.0 * std::max(ne_b[L], 0.0) * std::max(ne_b[R], 0.0)
                              / (std::max(ne_b[L], 0.0) + std::max(ne_b[R], 0.0) + 1e-300);
        fade_e = nmin / (nmin + drift_density);
      }
      const BoutReal heat =
          fade_e * (2.5 + ct) * ((current > 0.0) ? te_b[R] : te_b[L]) * current;
      const BoutReal work = fade_e * 0.5 * current * (phi_b[L] - phi_b[R]);
      // Current density along b at the face, for output
      const BoutReal jface = current / us.flux_over_B[f];
      const BoutReal w = std::abs(us.flux_over_B[f]);
      if (L < us.nlocal) {
        DivJpar(0, us.ys + L, 0) += current / us.vol[L];
        heating(0, us.ys + L, 0) += (heat + work) / us.vol[L];
        Jpar(0, us.ys + L, 0) += w * jface;
        weight(0, us.ys + L, 0) += w;
      }
      if (R < us.nlocal) {
        DivJpar(0, us.ys + R, 0) -= current / us.vol[R];
        heating(0, us.ys + R, 0) += (-heat + work) / us.vol[R];
        Jpar(0, us.ys + R, 0) += w * jface;
        weight(0, us.ys + R, 0) += w;
      }
    }
    for (int c = 0; c < us.nlocal; c++) {
      Jpar(0, us.ys + c, 0) /= std::max(weight(0, us.ys + c, 0), 1e-300);
    }
    if (electron_current) {
      add(electrons["energy_source"], heating);
    }
  }

  //////////////////////////////////////////////////////////////////
  // Drift advection of every charged species, and the drift current

  DivJdrift = zeroFrom(Ne);
  wall_particles = zeroFrom(Ne);
  wall_energy = zeroFrom(Ne);
  Sn_exb = zeroFrom(Ne);
  Sn_mag = zeroFrom(Ne);
  const std::vector<BoutReal> none(us.nfaces, 0.0);
  // Blend of second-order and upwind face values, from the electron density
  const std::vector<BoutReal> weight = us.secondOrderWeight(Ne);
  const std::vector<BoutReal> first_order(us.nfaces, 0.0);
  const std::vector<BoutReal>& weightE = second_order_exb ? weight : first_order;
  std::vector<BoutReal> flowE = exb_advection ? us.flowExB(phi) : none;
  // Drifts fade out in near-vacuum plasma: every drift flow through a face is
  // multiplied by n_f / (n_f + drift_density), n_f the harmonic mean of the
  // electron densities either side (the cell's own on a boundary face). The
  // GOAT grid fills the whole vessel, and the fast drift dynamics of its
  // near-vacuum far-SOL cells (n ~ 1e13 - 1e15 m^-3) stall the implicit
  // solver while carrying nothing that matters.
  std::vector<BoutReal> fade(us.nfaces, 1.0);
  if (drift_density > 0.0) {
    const auto ne_b = us.gather(Ne);
    for (int f = 0; f < us.nfaces; f++) {
      const BoutReal nL = std::max(ne_b[us.faceL[f]], 0.0);
      const BoutReal nR = std::max(ne_b[us.faceR[f]], 0.0);
      const BoutReal nf = (us.face_kind[f] == Unstructured::interior)
                              ? 2.0 * nL * nR / (nL + nR + 1e-300)
                              : nL;
      fade[f] = nf / (nf + drift_density);
      flowE[f] *= fade[f];
    }
  }
  if (exb_wall_density > 0.0) {
    // E x B flow into the wall fades out as the density of the cell falls
    // below exb_wall_density: otherwise it drains tenuous far-SOL cells
    // towards vacuum, where the implicit solver stalls. The flux it removes,
    // n v_E with n below about 1e16 m^-3, carries a negligible amount of
    // particles and energy.
    for (int f = 0; f < us.nfaces; f++) {
      if (us.face_kind[f] == Unstructured::wall) {
        const BoutReal nc = std::max(Ne(0, us.ys + us.faceL[f], 0), 0.0);
        flowE[f] *= nc / (nc + exb_wall_density);
      }
    }
  }

  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first];
    if (!species.isSet("charge") or (get<BoutReal>(species["charge"]) == 0.0)) {
      continue;
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    const BoutReal A = get<BoutReal>(species["AA"]);
    const Field3D N = getNoBoundary<Field3D>(species["density"]);
    const Field3D T = getNoBoundary<Field3D>(species["temperature"]);
    const Field3D P = species.isSet("pressure")
                          ? getNoBoundary<Field3D>(species["pressure"])
                          : Field3D{N * T};

    std::vector<BoutReal> flowD = magnetic_drift ? us.flowMagneticDrift(T, Z) : none;
    for (int f = 0; f < us.nfaces; f++) {
      flowD[f] *= fade[f];
    }
    std::vector<BoutReal> flowPE(us.nfaces), flowPD(us.nfaces);
    for (int f = 0; f < us.nfaces; f++) {
      flowPE[f] = 1.5 * flowE[f];
      flowPD[f] = 2.5 * flowD[f];
    }
    // Each drift carries the face value from its own upwind side. The magnetic
    // drift always gets the (blended) second-order values; E x B only with
    // second_order_exb
    const Field3D divD_N =
        magnetic_drift ? us.Div_f_flow(N, flowD, nullptr, &weight) : zeroFrom(N);
    if (magnetic_drift) {
      DivJdrift += Z * divD_N;
    }
    if (kv.first != "e") {
      // The electron density follows from quasineutrality.
      // Ions that the E x B drift carries into the wall: the flow through each
      // wall face goes into its ghost cell, for recycling and sputtering
      Field3D wall_flow{zeroFrom(N)};
      subtract(species["density_source"],
               us.Div_f_flow(N, flowE, &wall_flow, &weightE) + divD_N);
      set(species["wall_drift_flow"], wall_flow);
      if (Z > 0.0) {
        for (int f = 0; f < us.nfaces; f++) {
          if (us.face_kind[f] == Unstructured::wall) {
            wall_particles(0, us.ys + us.faceL[f], 0) +=
                wall_flow(0, us.ys + us.faceR[f], 0) / us.vol[us.faceL[f]];
          }
        }
      }
      if (species.isSet("velocity")) {
        const Field3D NV = N * getNoBoundary<Field3D>(species["velocity"]);
        subtract(species["momentum_source"],
                 A * (us.Div_f_flow(NV, flowE, nullptr, &weightE)
                      + us.Div_f_flow(NV, flowD, nullptr, &weight)));
      }
    }
    Field3D wall_heat{zeroFrom(N)}; // Energy out through each boundary face (wall and core)
    subtract(species["energy_source"], us.Div_f_flow(P, flowPE, &wall_heat, &weightE)
                                           + us.Div_f_flow(P, flowPD, &wall_heat, &weight));
    // For evolve_pressure's core_power_total, which reads the core ghosts only
    add(species["core_energy_flow"], wall_heat);
    for (int f = 0; f < us.nfaces; f++) {
      if (us.face_kind[f] == Unstructured::wall) {
        wall_energy(0, us.ys + us.faceL[f], 0) +=
            wall_heat(0, us.ys + us.faceR[f], 0) / us.vol[us.faceL[f]];
      }
    }
  }
  if (diagnose) {
    // E x B and magnetic drift particle sources of the electrons (= all ions)
    Sn_exb = -us.Div_f_flow(Ne, flowE, nullptr, &weightE);
    std::vector<BoutReal> flowDe = magnetic_drift ? us.flowMagneticDrift(Te, -1.0) : none;
    for (int f = 0; f < us.nfaces; f++) {
      flowDe[f] *= fade[f];
    }
    Sn_mag = magnetic_drift ? Field3D{-us.Div_f_flow(Ne, flowDe, nullptr, &weight)} : zeroFrom(Ne);
  }

  set(state["fields"]["phi_unstructured"], phi);
}

void DriftsUnstructured::finally(const Options& state) {
  const auto& us = Unstructured::get();
  const Options& allspecies = state["species"];
  const Options& electrons = allspecies["e"];
  const Field3D Ne = get<Field3D>(electrons["density"]);

  Field3D sigma_perp{zeroFrom(Ne)}, sigma_perp_over_n{zeroFrom(Ne)};
  for (int c = 0; c < us.nlocal; c++) {
    const int y = us.ys + c;
    sigma_perp(0, y, 0) = perp_mobility * softFloor(Ne(0, y, 0), 1e-7);
    sigma_perp_over_n(0, y, 0) = perp_mobility;
  }

  // Current into the sheath through wall faces: the ghost cells hold the
  // density and velocity of every species at the sheath edge
  DivJsheath = zeroFrom(Ne);
  for (const auto& kv : allspecies.getChildren()) {
    const Options& species = kv.second;
    if (!species.isSet("charge") or (get<BoutReal>(species["charge"]) == 0.0)
        or !species.isSet("velocity")) {
      continue;
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    const Field3D N = get<Field3D>(species["density"]);
    const Field3D V = get<Field3D>(species["velocity"]);
    for (int f = 0; f < us.nfaces; f++) {
      if (us.face_kind[f] != Unstructured::wall) {
        continue;
      }
      const int c = us.faceL[f];
      const int g = us.ys + us.faceR[f];
      DivJsheath(0, us.ys + c, 0) +=
          Z * N(0, g, 0) * V(0, g, 0) * us.flux_over_B[f] / us.vol[c];
    }
  }

  // Anomalous perpendicular current -sigma_perp (grad phi + grad p_i / n)
  DivJperp = -us.Div_a_Grad_perp(sigma_perp, phi, Unstructured::Scale::signed_);
  if (ion_pressure_term) {
    DivJperp -= us.Div_a_Grad_perp(sigma_perp_over_n, Pi_total);
  }

  ddt(phi) = -lambda * (DivJpar + DivJsheath + DivJdrift + DivJperp);
  us.zeroGhosts(ddt(phi));

  // For output, with the boundary values
  Ve = get<Field3D>(electrons["velocity"]);
  NVe = get<Field3D>(electrons["momentum"]);
}

void DriftsUnstructured::outputVars(Options& state) {
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  const BoutReal Jnorm = SI::qe * Nnorm * Cs0;
  const BoutReal rho_s0 = get<BoutReal>(state["rho_s0"]);

  state["phi"].setAttributes({{"time_dimension", "t"},
                              {"units", "V"},
                              {"conversion", Tnorm},
                              {"standard_name", "potential"},
                              {"long_name", "electrostatic potential relative to the wall"},
                              {"source", "drifts_unstructured"}});
  set_with_attrs(state["Ve"], Ve,
                 {{"time_dimension", "t"},
                  {"units", "m / s"},
                  {"conversion", Cs0},
                  {"long_name", "e parallel velocity"},
                  {"standard_name", "velocity"},
                  {"species", "e"},
                  {"source", "drifts_unstructured"}});
  set_with_attrs(state["Jpar"], Jpar,
                 {{"time_dimension", "t"},
                  {"units", "A m^-2"},
                  {"conversion", Jnorm},
                  {"long_name", "parallel current density"},
                  {"source", "drifts_unstructured"}});
  if (!diagnose) {
    return;
  }
  const std::vector<std::pair<std::string, Field3D*>> divs = {
      {"DivJpar", &DivJpar},
      {"DivJsheath", &DivJsheath},
      {"DivJdrift", &DivJdrift},
      {"DivJperp", &DivJperp}};
  for (const auto& [label, field] : divs) {
    set_with_attrs(state[label], *field,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", Jnorm / rho_s0},
                    {"long_name", "divergence of the current: " + label},
                    {"source", "drifts_unstructured"}});
  }
  set_with_attrs(state["Si_exb_wall"], wall_particles,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"long_name", "ions carried into the wall by the E x B drift, all ions"},
                  {"source", "drifts_unstructured"}});
  set_with_attrs(state["E_exb_wall"], wall_energy,
                 {{"time_dimension", "t"},
                  {"units", "W m^-3"},
                  {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                  {"long_name", "thermal energy carried into the wall by the E x B "
                                "drift, all species"},
                  {"source", "drifts_unstructured"}});
  set_with_attrs(state["Se_exb"], Sn_exb,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"long_name", "electron particle source of the E x B drift"},
                  {"source", "drifts_unstructured"}});
  set_with_attrs(state["Se_magdrift"], Sn_mag,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"long_name", "electron particle source of the magnetic drift"},
                  {"source", "drifts_unstructured"}});
}
