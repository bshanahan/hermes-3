#pragma once
#ifndef WALL_SPUTTERING_H
#define WALL_SPUTTERING_H

#include "component.hxx"
#include "unstructured.hxx"

#include <bout/constants.hxx>
#include <bout/utils.hxx>

#include <cmath>
#include <fmt/format.h>
#include <string>
#include <vector>

/// Sputtering of wall material with a fixed yield, on a flat-cell
/// (unstructured) grid.
///
/// Every wall face releases `yield` atoms of the sputtered species for each
/// particle of the impinging species that reaches it:
///   ions      the flow through the sheath (ghost cell values), plus the
///             cross-field loss if anomalous_diffusion has a wall_decay_length
///   neutrals  the one-way thermal flux n * sqrt(8 T / (pi m)) / 4 on the face
/// The atoms enter the cell next to the face with a fixed energy each.
/// The yield does not depend on energy, angle or surface temperature, as for
/// the chemical sputtering yield of SOLPS-ITER carbon cases.
///
/// Options, in the section named after the component:
///   species   Sputtered neutral species, e.g. c
///   from      Impinging species, e.g. "d+, d"
///   yield     Atoms released per impinging particle
///   energy    Energy of each released atom [eV] (default 1 eV)
///   boxes     Optional: only wall faces in these boxes sputter,
///             "Rmin, Rmax, Zmin, Zmax | ..." in metres
struct WallSputtering : public NamedComponent<WallSputtering> {
  WallSputtering(std::string name, Options& alloptions, Solver*)
      : NamedComponent(name, {readOnly("species:{from}:{inputs}"),
                              readIfSet("species:{from}:{optional}"),
                              readWrite("species:{sp}:{outputs}")}) {
    if (!Unstructured::active()) {
      throw BoutException("wall_sputtering needs a flat-cell (unstructured) grid");
    }
    const auto& us = Unstructured::get();
    Options& options = alloptions[name];
    const auto& units = alloptions["units"];
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);
    const BoutReal Tnorm = get<BoutReal>(units["eV"]);

    species = trim(options["species"].doc("Sputtered neutral species").as<std::string>());
    for (const auto& item :
         strsplit(options["from"].doc("Impinging species, e.g. 'd+, d'").as<std::string>(),
                  ',')) {
      const std::string sp = trim(item, " \t()");
      if (!sp.empty()) {
        from.push_back(sp);
      }
    }
    yield = options["yield"]
                .doc("Atoms released per impinging particle")
                .withDefault<BoutReal>(0.0);
    energy = options["energy"]
                 .doc("Energy of each released atom [eV]")
                 .withDefault<BoutReal>(1.0)
             / Tnorm;
    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);

    active_face.assign(us.nfaces, false);
    if (options.isSet("boxes")) {
      for (const int f : us.wallFacesInBoxes(options["boxes"].as<std::string>())) {
        active_face[f] = true;
      }
    } else {
      for (int f = 0; f < us.nfaces; f++) {
        active_face[f] = us.face_kind[f] == Unstructured::wall;
      }
    }

    substitutePermissions("sp", {species});
    substitutePermissions("from", from);
    substitutePermissions("inputs", {"AA", "density", "temperature"});
    substitutePermissions("optional", {"velocity", "charge", "wall_loss_flow"});
    substitutePermissions("outputs", {"density_source", "energy_source"});
  };

  void outputVars(Options& state) override {
    if (diagnose) {
      set_with_attrs(state[fmt::format("S{}_sputter", species)], particle_source,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"long_name", "sputtered particle source"},
                      {"source", "wall_sputtering"}});
    }
  }

  static constexpr auto type = "wall_sputtering";

private:
  std::string species;
  std::vector<std::string> from;
  BoutReal Nnorm, Omega_ci;
  BoutReal yield;
  BoutReal energy; ///< Energy per released atom, normalised
  bool diagnose;
  std::vector<bool> active_face;
  Field3D particle_source; ///< Normalised particle source

  void transform_impl(GuardedOptions& state) override {
    const auto& us = Unstructured::get();
    particle_source = 0.0;
    particle_source.allocate();
    for (const auto& name : from) {
      auto sp = state["species"][name];
      const Field3D N = get<Field3D>(sp["density"]);
      const Field3D T = get<Field3D>(sp["temperature"]);
      const BoutReal AA = get<BoutReal>(sp["AA"]);
      const bool ion = sp.isSet("charge") and (get<BoutReal>(sp["charge"]) != 0.0);
      const bool has_velocity = sp.isSet("velocity");
      const Field3D V = has_velocity ? get<Field3D>(sp["velocity"]) : Field3D{0.0};
      const bool has_loss = ion and sp.isSet("wall_loss_flow");
      const Field3D loss = has_loss ? get<Field3D>(sp["wall_loss_flow"]) : Field3D{0.0};
      for (int f = 0; f < us.nfaces; f++) {
        if (!active_face[f]) {
          continue;
        }
        const int c = us.faceL[f];
        const int i = us.ys + c;
        const int g = us.ys + us.faceR[f];
        BoutReal flow = 0.0;
        if (ion) {
          // Flow out through the sheath, from the ghost cell values
          if (has_velocity) {
            flow = std::max(N(0, g, 0) * V(0, g, 0) * us.flux_over_B[f], 0.0);
          }
          if (has_loss) {
            flow += loss(0, g, 0);
          }
        } else {
          const BoutReal tn = std::max(T(0, i, 0), 0.0);
          flow = 0.25 * sqrt(8 * tn / (PI * AA)) * std::max(N(0, i, 0), 0.0)
                 * us.face_area[f];
        }
        particle_source(0, i, 0) += yield * flow / us.vol[c];
      }
    }
    auto to = state["species"][species];
    add(to["density_source"], particle_source);
    add(to["energy_source"], energy * particle_source);
  };
};

namespace {
RegisterComponent<WallSputtering> register_wall_sputtering;
}

#endif // WALL_SPUTTERING_H
