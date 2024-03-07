#pragma once
#ifndef EXTERNAL_FIELD_H
#define EXTERNAL_FIELD_H

#include "component.hxx"

/// Calculate diamagnetic flows

struct ExternalField : public Component {
  ExternalField(std::string name, Options &options, Solver *UNUSED(solver));

  /// For every species, if it has:
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  BoutReal Lnorm;

  Field3D psiext;
  Field3D beta_em;
  bool diagnose;
  bool psiext_time_dependent;
  FieldGeneratorPtr psiext_factory;
  bool time_normalisation;
  std::string source_string;

  Options& options;

};

namespace {
RegisterComponent<ExternalField> registercomponentexternalfield("external_field");
}

#endif // EXTERNAL_FIELD_H


