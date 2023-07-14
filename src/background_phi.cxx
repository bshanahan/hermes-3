
#include "../include/background_phi.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

BackgroundPhi::BackgroundPhi(std::string name, Options& alloptions, Solver* solver)
  : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  phi_bg = options["phi_bg"]
    .doc("Background potential.")
    .withDefault(Field3D(0.));
  
}

void BackgroundPhi::transform(Options& state) {
  AUTO_TRACE();
  Field3D phi = get<Field3D>(state["fields"]["phi"]);

  // add the background due to the imposed radial electric field.
  set(state["fields"]["phi"],phi + phi_bg);
  
}



