#include "VEMassTimeDerivative.h"

registerMooseObject("wombatApp", VEMassTimeDerivative);

InputParameters
VEMassTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addClassDescription(
      "Depth-integrated mass storage term for one fluid phase in a vertical-equilibrium "
      "simulation: d/dt(H * phi_bar * rho_c * S_c_bar). Jacobian computed automatically via AD.");

  params.addParam<unsigned int>(
      "fluid_phase",
      0,
      "Fluid phase index this kernel acts on. 0 = CO2 (non-wetting), 1 = brine (wetting). "
      "Use fluid_phase=0 on variable=sat_n and fluid_phase=1 on variable=pp_top.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEMassTimeDerivative::VEMassTimeDerivative(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _H(getMaterialProperty<Real>("ve_H")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _density_old(getMaterialPropertyOld<std::vector<Real>>("ve_density")),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _saturation_old(getMaterialPropertyOld<std::vector<Real>>("ve_saturation"))
{
  // VE is always a two-phase CO2-brine system (phase 0 = CO2, 1 = brine).
  if (_fluid_phase > 1)
    paramError("fluid_phase",
               "fluid_phase=",
               _fluid_phase,
               " but a VE simulation has only 2 fluid phases (0 = CO2, 1 = brine).");
}

ADReal
VEMassTimeDerivative::precomputeQpResidual()
{
  // New accumulation: AD type so density derivatives (wrt pp_top) and
  // saturation derivatives (wrt sat_n) are tracked automatically.
  const ADReal accumulation_new =
      _H[_qp] * _phi_bar[_qp] * _density[_qp][_fluid_phase] * _saturation[_qp][_fluid_phase];

  // Old accumulation: plain Real -- no AD derivatives, no Jacobian contribution.
  const Real accumulation_old = _H[_qp] * _phi_bar[_qp] * _density_old[_qp][_fluid_phase] *
                                _saturation_old[_qp][_fluid_phase];

  return (accumulation_new - accumulation_old) / _dt;
}
