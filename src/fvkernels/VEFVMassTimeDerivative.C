#include "VEFVMassTimeDerivative.h"

registerMooseObject("wombatApp", VEFVMassTimeDerivative);

InputParameters
VEFVMassTimeDerivative::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Finite-volume depth-integrated mass storage term for one fluid phase: "
      "d/dt(H * phi_bar * rho_c * S_c). FV analogue of VEMassTimeDerivative.");

  params.addRequiredParam<UserObjectName>("VEDictator", "The VEDictator UserObject.");
  params.addParam<unsigned int>("fluid_phase", 0,
      "Fluid phase index (0 = CO2, 1 = brine). Use fluid_phase=0 on variable=sat_n "
      "and fluid_phase=1 on variable=pp_top.");
  params.addRequiredCoupledVar("z_top",
      "Top-surface elevation z_T [m] as an FV variable. H = z_top - z_bottom.");
  params.addRequiredCoupledVar("z_bottom",
      "Bottom-surface elevation z_B [m] as an FV variable. H = z_top - z_bottom.");
  return params;
}

VEFVMassTimeDerivative::VEFVMassTimeDerivative(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _dictator(getUserObject<VEDictator>("VEDictator")),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _z_top(coupledValue("z_top")),
    _z_bottom(coupledValue("z_bottom")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _density_old(getMaterialPropertyOld<std::vector<Real>>("ve_density")),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _saturation_old(getMaterialPropertyOld<std::vector<Real>>("ve_saturation"))
{
  if (_fluid_phase >= _dictator.numPhases())
    paramError("fluid_phase", "fluid_phase=", _fluid_phase,
               " but VEDictator reports only ", _dictator.numPhases(), " phases.");
}

ADReal
VEFVMassTimeDerivative::computeQpResidual()
{
  const Real H = _z_top[_qp] - _z_bottom[_qp];

  const ADReal accum_new = H * _phi_bar[_qp]
                         * _density[_qp][_fluid_phase]
                         * _saturation[_qp][_fluid_phase];

  const Real accum_old = H * _phi_bar[_qp]
                       * _density_old[_qp][_fluid_phase]
                       * _saturation_old[_qp][_fluid_phase];

  return (accum_new - accum_old) / _dt;
}
