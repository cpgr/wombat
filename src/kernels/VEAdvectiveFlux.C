#include "VEAdvectiveFlux.h"

// Not registered -- VEAdvectiveFluxBase is abstract.
// Use VEAdvectiveFluxP or VEAdvectiveFluxS in input files.

InputParameters
VEAdvectiveFluxBase::validParams()
{
  InputParameters params = ADKernelGrad::validParams();

  params.addParam<unsigned int>(
      "fluid_phase",
      0,
      "Fluid phase index. 0 = CO2 (non-wetting), 1 = brine (wetting).");

  RealVectorValue g_default(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>(
      "gravity",
      g_default,
      "Gravity vector (m/s2). Defaults to (0, 0, -9.81). "
      "Only the magnitude |g| enters the depth-integrated buoyancy term.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEAdvectiveFluxBase::VEAdvectiveFluxBase(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _H(getMaterialProperty<Real>("ve_H")),
    _K_up(getMaterialProperty<RealTensorValue>("ve_K_up")),
    _relperm(getADMaterialProperty<std::vector<Real>>("ve_relperm")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("ve_viscosity")),
    _grad_z_top(getMaterialProperty<RealGradient>("ve_grad_z_top"))
{
  // VE is always a two-phase CO2-brine system (phase 0 = CO2, 1 = brine).
  if (_fluid_phase > 1)
    paramError("fluid_phase",
               "fluid_phase=",
               _fluid_phase,
               " but a VE simulation has only 2 fluid phases (0 = CO2, 1 = brine).");

  if (_gravity_magnitude == 0.0)
    paramWarning("gravity",
                 "Gravity magnitude is zero. The buoyancy term rho_c * |g| * grad(z_T) "
                 "will be absent. This is correct only for flat-formation verification cases.");
}
