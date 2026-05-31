#include "VEUpscaledCapPressure.h"

registerMooseObject("wombatApp", VEUpscaledCapPressure);

InputParameters
VEUpscaledCapPressure::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Upscaled column capillary pressure ve_pc_up = (rho_w - rho_n) * |g| * h + pc_entry, "
      "with h the reconstructed plume thickness (ve_h). Diagnostic (non-AD) for now, "
      "paralleling ve_h.");

  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only the magnitude enters the buoyancy term.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Capillary entry/fringe pressure added at the plume top [Pa].");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEUpscaledCapPressure::VEUpscaledCapPressure(const InputParameters & parameters)
  : Material(parameters),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _h(getMaterialProperty<Real>("ve_h")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _pc_up(declareProperty<Real>("ve_pc_up"))
{
}

void
VEUpscaledCapPressure::computeQpProperties()
{
  // delta_rho = rho_w - rho_n (brine minus CO2). Density carries no Jacobian
  // here because ve_pc_up is a diagnostic Real; take raw values.
  const Real delta_rho =
      MetaPhysicL::raw_value(_density[_qp][1]) - MetaPhysicL::raw_value(_density[_qp][0]);

  _pc_up[_qp] = delta_rho * _gravity_magnitude * _h[_qp] + _pc_entry;
}
