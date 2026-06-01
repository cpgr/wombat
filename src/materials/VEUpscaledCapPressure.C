#include "VEUpscaledCapPressure.h"

registerMooseObject("wombatApp", VEUpscaledCapPressure);

InputParameters
VEUpscaledCapPressure::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Upscaled column capillary pressure ve_pc_up = (rho_w - rho_n) * |g| * h + pc_entry "
      "and its partial derivatives ve_dPcup_dsatn = (rho_w-rho_n)*|g|*H/(1-S_wr) and "
      "ve_dPcup_dH = (rho_w-rho_n)*|g|*sat_n/(1-S_wr). "
      "ve_pc_up and ve_dPcup_dH are AD properties (Jacobian flows through density / sat_n). "
      "Both derivatives are consumed by VEAdvectiveFluxS when capillary=true to assemble "
      "grad(Pc^up) = ve_dPcup_dsatn*grad(sat_n) + ve_dPcup_dH*grad(H).");

  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only the magnitude enters the buoyancy term.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Capillary entry/fringe pressure added at the plume top [Pa].");
  params.addRangeCheckedParam<Real>(
      "S_wr",
      0.0,
      "S_wr >= 0 & S_wr < 1",
      "Residual water saturation [-]. Must match the value used in VEPlumeReconstruction "
      "so that ve_dPcup_dsatn = (rho_w-rho_n)*|g|*H/(1-S_wr) is consistent with ve_h.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEUpscaledCapPressure::VEUpscaledCapPressure(const InputParameters & parameters)
  : Material(parameters),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _S_wr(getParam<Real>("S_wr")),
    _h(getADMaterialProperty<Real>("ve_h")),
    _density(getADMaterialProperty<std::vector<Real>>("ve_density")),
    _H(getMaterialProperty<Real>("ve_H")),
    _pc_up(declareADProperty<Real>("ve_pc_up")),
    _dPcup_dsatn(declareProperty<Real>("ve_dPcup_dsatn")),
    _dPcup_dH(declareADProperty<Real>("ve_dPcup_dH"))
{
}

void
VEUpscaledCapPressure::computeQpProperties()
{
  const ADReal delta_rho = _density[_qp][1] - _density[_qp][0];
  _pc_up[_qp] = delta_rho * _gravity_magnitude * _h[_qp] + _pc_entry;

  // Sharp-interface formula: d(Pc^up)/d(sat_n) = delta_rho * g * H / (1 - S_wr).
  // Raw value used because density is constant in current fluid property models;
  // the grad(sat_n) term's Jacobian flows through the variable gradient _grad_u.
  const Real delta_rho_val = MetaPhysicL::raw_value(delta_rho);
  _dPcup_dsatn[_qp] = delta_rho_val * _gravity_magnitude * _H[_qp] / (1.0 - _S_wr);

  // Sharp-interface formula: d(Pc^up)/d(H) = delta_rho * g * sat_n / (1 - S_wr)
  // = delta_rho * g * h / H  (since sharp h = sat_n*H/(1-S_wr)). This multiplies
  // grad(H) in VEAdvectiveFluxS, completing grad(Pc^up) = delta_rho*g*grad(h)
  // for laterally varying thickness. It is AD because grad(H) is fixed geometry
  // (no variable gradient to carry the Jacobian), so the d/d(sat_n) derivative
  // must ride on this coefficient: d(ve_dPcup_dH)/d(sat_n) = delta_rho*g/(1-S_wr).
  _dPcup_dH[_qp] = delta_rho * _gravity_magnitude * _h[_qp] / _H[_qp];
}
