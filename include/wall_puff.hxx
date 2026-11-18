#pragma once
#ifndef WALL_PUFF_H
#define WALL_PUFF_H

#include "component.hxx"
#include "unstructured.hxx"

#include <bout/boutcomm.hxx>
#include <bout/constants.hxx>
#include <bout/utils.hxx>

#include <cmath>
#include <fmt/format.h>
#include <string>
#include <vector>

/// Gas puff from valves on the wall of a flat-cell (unstructured) grid.
///
/// Each valve is given by a position (R, Z) in metres and is moved to the
/// nearest wall face. Its gas enters the cell next to that face as a source
/// of neutral particles, with a fixed energy per particle and, optionally,
/// momentum along the wall normal (only the component along the magnetic
/// field can be carried by the neutral model).
///
/// Options, in the section named after the component:
///   species    Neutral species to puff, e.g. d
///   locations  "R1, Z1 | R2, Z2 | ..."  [m]
///   rates      "r1 | r2 | ..."  particles per second from each valve, or
///   total_rate particles per second shared equally between the valves
///   energy     Energy per particle [eV] (default 3 eV)
struct WallPuff : public NamedComponent<WallPuff> {
  WallPuff(std::string name, Options& alloptions, Solver*)
      : NamedComponent(name, {readOnly("species:{sp}:AA"),
                              readWrite("species:{sp}:{outputs}")}) {
    if (!Unstructured::active()) {
      throw BoutException("wall_puff needs a flat-cell (unstructured) grid");
    }
    const auto& us = Unstructured::get();
    Options& options = alloptions[name];
    const auto& units = alloptions["units"];
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);
    const BoutReal Tnorm = get<BoutReal>(units["eV"]);

    species = trim(options["species"].doc("Neutral species to puff").as<std::string>());
    energy = options["energy"]
                 .doc("Energy of each puffed particle [eV]")
                 .withDefault<BoutReal>(3.0)
             / Tnorm;
    wall_normal_launch =
        options["wall_normal_launch"]
            .doc("Give the gas momentum along the wall normal? Only its component "
                 "along the magnetic field is used")
            .withDefault<bool>(true);
    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);

    // Valve positions
    std::vector<std::pair<BoutReal, BoutReal>> positions;
    for (const auto& item :
         strsplit(options["locations"]
                      .doc("Valve positions 'R, Z | R, Z | ...' in metres")
                      .as<std::string>(),
                  '|')) {
      if (trim(item).empty()) {
        continue;
      }
      const auto parts = strsplit(item, ',');
      const std::vector<std::string> rz(parts.begin(), parts.end());
      if (rz.size() != 2) {
        throw BoutException("wall_puff location '{}' should be 'R, Z' in metres", item);
      }
      positions.emplace_back(stringToReal(trim(rz[0])), stringToReal(trim(rz[1])));
    }
    if (positions.empty()) {
      throw BoutException("wall_puff '{}' has no locations", name);
    }

    // Rates [particles / s]
    std::vector<BoutReal> rates;
    if (options.isSet("rates")) {
      for (const auto& item : strsplit(
               options["rates"].doc("Particles per second from each valve").as<std::string>(),
               '|')) {
        if (!trim(item).empty()) {
          rates.push_back(stringToReal(trim(item)));
        }
      }
      if (rates.size() != positions.size()) {
        throw BoutException("wall_puff '{}': {} rates for {} locations", name,
                            rates.size(), positions.size());
      }
    } else {
      const BoutReal total =
          options["total_rate"]
              .doc("Particles per second, shared equally between the valves")
              .withDefault<BoutReal>(0.0);
      rates.assign(positions.size(), total / static_cast<BoutReal>(positions.size()));
    }

    // Move each valve to the nearest wall face, which may be on any processor
    int rank = 0;
    MPI_Comm_rank(BoutComm::get(), &rank);
    particle_source = 0.0;
    launch_factor = 0.0;
    particle_source.allocate();
    launch_factor.allocate();
    for (std::size_t k = 0; k < positions.size(); k++) {
      struct {
        double distance;
        int rank;
      } mine{1e300, rank}, best{};
      int face = -1;
      for (int f = 0; f < us.nfaces; f++) {
        if (us.face_kind[f] != Unstructured::wall) {
          continue;
        }
        const double d = std::hypot(us.face_R[f] - positions[k].first,
                                    us.face_Z[f] - positions[k].second);
        if (d < mine.distance) {
          mine.distance = d;
          face = f;
        }
      }
      MPI_Allreduce(&mine, &best, 1, MPI_DOUBLE_INT, MPI_MINLOC, BoutComm::get());
      if (best.rank != rank) {
        continue; // Another processor owns the nearest face
      }
      const int c = us.faceL[face];
      // Particles per second -> per radian, per volume, normalised
      const BoutReal source = rates[k] / (TWOPI * us.vol[c] * std::pow(us.length_norm, 3))
                              / (Nnorm * Omega_ci);
      particle_source(0, us.ys + c, 0) += source;
      // b . (normal into the plasma)
      launch_factor(0, us.ys + c, 0) -= source * us.flux_over_B[face] / us.face_area[face];
      output.write("\twall_puff {}: valve at (R, Z) = ({:.3f}, {:.3f}) m -> wall face at "
                   "({:.3f}, {:.3f}) m, {:.1f} mm away, {:e} particles/s\n",
                   name, positions[k].first, positions[k].second, us.face_R[face],
                   us.face_Z[face], 1e3 * best.distance, rates[k]);
    }

    substitutePermissions("sp", {species});
    substitutePermissions("outputs",
                          {"density_source", "energy_source", "momentum_source"});
  };

  void outputVars(Options& state) override {
    if (diagnose) {
      set_with_attrs(state[fmt::format("S{}_puff", species)], particle_source,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"long_name", "gas puff particle source"},
                      {"source", "wall_puff"}});
    }
  }

  static constexpr auto type = "wall_puff";

private:
  std::string species;
  BoutReal Nnorm, Omega_ci;
  BoutReal energy;         ///< Energy per particle, normalised
  bool wall_normal_launch; ///< Add momentum along the wall normal?
  bool diagnose;
  Field3D particle_source; ///< Normalised particle source
  Field3D launch_factor;   ///< Particle source times b . (inward wall normal)

  void transform_impl(GuardedOptions& state) override {
    auto sp = state["species"][species];
    add(sp["density_source"], particle_source);
    add(sp["energy_source"], energy * particle_source);
    if (wall_normal_launch) {
      const BoutReal AA = get<BoutReal>(sp["AA"]);
      // Mean wall-normal speed of a cosine distribution at this energy
      const BoutReal speed = (2. / 3) * sqrt(2. * energy / AA);
      add(sp["momentum_source"], (AA * speed) * launch_factor);
    }
  };
};

namespace {
RegisterComponent<WallPuff> register_wall_puff;
}

#endif // WALL_PUFF_H
