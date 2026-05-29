#include "VEAdvectiveFluxS.h"

registerMooseObject("wombatApp", VEAdvectiveFluxS);

InputParameters
VEAdvectiveFluxS::validParams()
{
  InputParameters params = VEAdvectiveFluxBase::validParams();
  params.addClassDescription(
      "Depth-integrated Darcy advective flux for the saturation variable (variable=sat_n). "
      "Couples to pp_top to obtain grad(pp_top) for Darcy flow. Intended for the CO2 mass "
      "equation (fluid_phase=0). Jacobian (including off-diagonal wrt pp_top) via AD.");

  params.addRequiredCoupledVar("pp_top",
                               "The pressure-at-top-surface primary variable. "
                               "Its gradient is used as the Darcy pressure gradient.");
  return params;
}

VEAdvectiveFluxS::VEAdvectiveFluxS(const InputParameters & parameters)
  : VEAdvectiveFluxBase(parameters),
    _grad_pp_top(adCoupledGradient("pp_top"))
{
}

ADRealVectorValue
VEAdvectiveFluxS::precomputeQpResidual()
{
  const ADReal & rho_c = _density[_qp][_fluid_phase];
  const ADReal & kr_c  = _relperm[_qp][_fluid_phase];
  const ADReal & mu_c  = _viscosity[_qp][_fluid_phase];

  // _grad_pp_top[_qp] is adCoupledGradient("pp_top"), carrying AD derivatives
  // wrt pp_top DOFs so the off-diagonal Jacobian block is assembled correctly.
  const ADRealVectorValue potential_gradient =
      _grad_pp_top[_qp] + rho_c * _gravity_magnitude * _grad_z_top[_qp];

  const ADReal mobility = _H[_qp] * kr_c * rho_c / mu_c;

  return mobility * (_K_up[_qp] * potential_gradient);
}
