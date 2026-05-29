#include "VEAdvectiveFluxP.h"

registerMooseObject("wombatApp", VEAdvectiveFluxP);

InputParameters
VEAdvectiveFluxP::validParams()
{
  InputParameters params = VEAdvectiveFluxBase::validParams();
  params.addClassDescription(
      "Depth-integrated Darcy advective flux for the pressure variable (variable=pp_top). "
      "Uses _grad_u as grad(pp_top) directly. Intended for the brine mass equation "
      "(fluid_phase=1). Jacobian assembled automatically via AD.");
  return params;
}

VEAdvectiveFluxP::VEAdvectiveFluxP(const InputParameters & parameters)
  : VEAdvectiveFluxBase(parameters)
{
}

ADRealVectorValue
VEAdvectiveFluxP::precomputeQpResidual()
{
  const ADReal & rho_c = _density[_qp][_fluid_phase];
  const ADReal & kr_c  = _relperm[_qp][_fluid_phase];
  const ADReal & mu_c  = _viscosity[_qp][_fluid_phase];

  // _grad_u[_qp] is grad(pp_top) because this kernel acts on pp_top.
  const ADRealVectorValue potential_gradient =
      _grad_u[_qp] + rho_c * _gravity_magnitude * _grad_z_top[_qp];

  const ADReal mobility = _H[_qp] * kr_c * rho_c / mu_c;

  return mobility * (_K_up[_qp] * potential_gradient);
}
