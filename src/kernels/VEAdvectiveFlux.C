#include "VEAdvectiveFlux.h"

// Not registered -- VEAdvectiveFluxBase is abstract.
// Use VEAdvectiveFluxP or VEAdvectiveFluxS in input files.

InputParameters
VEAdvectiveFluxBase::validParams()
{
  InputParameters params = ADKernelGrad::validParams();

  params.addRequiredParam<UserObjectName>("VEDictator",
                                         "The VEDictator UserObject for this simulation.");

  params.addParam<unsigned int>(
      "fluid_phase",
      0,
      "Fluid phase index. 0 = CO2 (non-wetting), 1 = brine (wetting).");

  RealVectorValue g_default(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>(
      "gravity",
      g_default,
      "Gravity vector [m/s2]. Defaults to (0, 0, -9.81). "
      "Only the magnitude |g| enters the depth-integrated buoyancy term.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEAdvectiveFluxBase::VEAdvectiveFluxBase(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _dictator(getUserObject<VEDictator>("VEDictator")),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _H(getMaterialProperty<Real>("ve_H")),
    _K_up(getMaterialProperty<RealTensorValue>("ve_K_up")),
    _relperm(getADMaterialProperty<std::vector<Real>>("ve_relperm")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("ve_viscosity")),
    _grad_z_top(getMaterialProperty<RealGradient>("ve_grad_z_top"))
{
  if (_fluid_phase >= _dictator.numPhases())
    paramError("fluid_phase",
               "fluid_phase=",
               _fluid_phase,
               " but the VEDictator reports only ",
               _dictator.numPhases(),
               " fluid phases (0-indexed).");

  if (_gravity_magnitude == 0.0)
    paramWarning("gravity",
                 "Gravity magnitude is zero. The buoyancy term rho_c * |g| * grad(z_T) "
                 "will be absent. This is correct only for flat-formation verification cases.");
}
