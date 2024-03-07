#include "../include/external_field.hxx"
#include <bout/constants.hxx>
#include <bout/difops.hxx>
#include <bout/field_factory.hxx>

using bout::globals::mesh;

ExternalField::ExternalField(std::string name, Options& alloptions,
			     Solver* UNUSED(solver)): options{alloptions[name]} {

  psiext_time_dependent = options["time_dependent"]
    .doc("Use a time-dependent external field?")
    .withDefault<bool>(true);
  
  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  Lnorm = get<BoutReal>(units["meters"]);
  BoutReal Bnorm = units["Tesla"];
  BoutReal Tnorm = units["eV"];
  BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  beta_em = SI::mu0 * SI::qe * Tnorm * Nnorm / SQ(Bnorm);

  time_normalisation = 1./Omega_ci;
  source_string = options["function"]
      .doc("Time-dependent external field.")
      .as<std::string>();
  if (psiext_time_dependent) {
    psiext_factory = FieldFactory::get()->parse(source_string, &options);
  } else {
    psiext =
      options["function"].doc("Form of the external field").withDefault<Field3D>(0.0);
    psiext /= Lnorm;
  }
}

void ExternalField::transform(Options& state) {
  // Iterate through all subsections
  Options& allspecies = state["species"];

  if (psiext_time_dependent) {
    // Evaluate the source_prefactor function at the current time in seconds and scale source with it
    BoutReal time = get<BoutReal>(state["time"]) * time_normalisation;
    psiext = FieldFactory::get()->create3D(source_string, &options, nullptr, CELL_CENTRE, time); //psiext_factory ->generate(bout::generator::Context().set("t",time*time_normalisation));
    psiext /= Lnorm;
  }
  
  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }

    if (IS_SET(state["fields"]["Apar"])){
      const Field3D Apar = state["fields"]["Apar"];
      psiext += Apar;
    }
    
    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);
      subtract(species["density_source"], 0.5*beta_em*bracket(psiext, N, BRACKET_ARAKAWA));
    }

    if (IS_SET(species["vorticity"])) {
      auto Vort = GET_VALUE(Field3D, species["vorticity"]);
      subtract(species["vorticity_source"], 0.5*beta_em*bracket(psiext, Vort, BRACKET_ARAKAWA));
    }

    if (IS_SET(species["momentum"])) {
      auto NV = GET_VALUE(Field3D, species["momentum"]);
      subtract(species["momentum_source"], 0.5*beta_em*bracket(psiext, NV, BRACKET_ARAKAWA));
    }
     
    if (IS_SET(species["pressure"])) {
      auto P = GET_VALUE(Field3D, species["pressure"]);
      subtract(species["energy_source"], 0.5*beta_em*bracket(psiext, P, BRACKET_ARAKAWA));
    }

  }
}

void ExternalField::outputVars(Options &state) {
  // Normalisations
  auto Bnorm = get<BoutReal>(state["Bnorm"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  set_with_attrs(state["beta_em"], beta_em, {
      {"long_name", "Helmholtz equation parameter"}
    });

  set_with_attrs(state["psiext"], psiext, {
      {"time_dimension", "t"},
      {"units", "T m"},
      {"conversion", Bnorm * rho_s0},
      {"long_name", "External, applied field."}
    });
}
