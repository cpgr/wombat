#include "VEAdvectiveFluxS.h"

registerMooseObject("wombatApp", VEAdvectiveFluxS);

InputParameters
VEAdvectiveFluxS::validParams()
{
  InputParameters params = VEAdvectiveFluxBase::validParams();
  params.addClassDescription(
      "Depth-integrated Darcy advective flux for the saturation variable (variable=sat_n). "
      "Couples to pp_top to obtain grad(pp_top) for Darcy flow. Intended for the CO2 mass "
      "equation (fluid_phase=0). Jacobian (including off-diagonal wrt pp_top) via AD. "
      "With capillary=true, adds the capillary-diffusion term grad(Pc^up) = "
      "ve_dPcup_dsatn * grad(sat_n) + ve_dPcup_dH * grad(H) to the Darcy potential; "
      "requires VEUpscaledCapPressure (and VEGeometry's ve_grad_H).");

  params.addRequiredCoupledVar("pp_top",
                               "The pressure-at-top-surface primary variable. "
                               "Its gradient is used as the Darcy pressure gradient.");

  params.addParam<bool>(
      "capillary",
      false,
      "If true, adds grad(Pc^up) = d(Pc^up)/d(sat_n)*grad(sat_n) + "
      "d(Pc^up)/d(H)*grad(H) to the CO2 Darcy potential (two-pressure VE). "
      "Requires VEUpscaledCapPressure (declares ve_dPcup_dsatn, ve_dPcup_dH) and "
      "VEGeometry (declares ve_grad_H). The grad(H) term vanishes for "
      "constant-thickness formations, so existing flat/constant-H golds are "
      "unchanged. Default OFF so existing sharp-interface inputs are unchanged.");

  return params;
}

VEAdvectiveFluxS::VEAdvectiveFluxS(const InputParameters & parameters)
  : VEAdvectiveFluxBase(parameters),
    _grad_pp_top(adCoupledGradient("pp_top")),
    _capillary(getParam<bool>("capillary")),
    _dPcup_dsatn(_capillary ? &getMaterialProperty<Real>("ve_dPcup_dsatn") : nullptr),
    _dPcup_dH(_capillary ? &getADMaterialProperty<Real>("ve_dPcup_dH") : nullptr),
    _grad_H(_capillary ? &getMaterialProperty<RealGradient>("ve_grad_H") : nullptr)
{
}

ADRealVectorValue
VEAdvectiveFluxS::precomputeQpResidual()
{
  const ADReal & rho_c = _density[_qp][_fluid_phase];
  const ADReal & kr_c  = _relperm[_qp][_fluid_phase];
  const ADReal & mu_c  = _viscosity[_qp][_fluid_phase];

  // _grad_pp_top[_qp] carries AD derivatives wrt pp_top DOFs (off-diagonal Jacobian).
  // _grad_u[_qp]      carries AD derivatives wrt sat_n DOFs  (diagonal block).
  ADRealVectorValue potential_gradient =
      _grad_pp_top[_qp] + rho_c * _gravity_magnitude * _grad_z_top[_qp];

  if (_capillary)
  {
    // grad(Pc^up) = d(Pc^up)/d(sat_n) * grad(sat_n) + d(Pc^up)/d(H) * grad(H).
    // The first term's sat_n Jacobian rides on the AD variable gradient _grad_u;
    // the second uses the fixed-geometry ve_grad_H, so its sat_n Jacobian must
    // come from the AD coefficient _dPcup_dH. Together they reconstruct the full
    // delta_rho*g*grad(h) capillary drive for laterally varying thickness.
    potential_gradient += (*_dPcup_dsatn)[_qp] * _grad_u[_qp];
    potential_gradient += (*_dPcup_dH)[_qp] * (*_grad_H)[_qp];
  }

  const ADReal mobility = _H[_qp] * kr_c * rho_c / mu_c;

  return mobility * (_K_up[_qp] * potential_gradient);
}
