#include "VEUpscaledCapPressure.h"

registerMooseObject("wombatApp", VEUpscaledCapPressure);

InputParameters
VEUpscaledCapPressure::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Upscaled column capillary pressure ve_pc_up = (rho_w - rho_n) * |g| * h + pc_entry "
      "and its saturation derivative ve_dPcup_dsatn = (rho_w-rho_n)*|g|*H/(1-S_wr). "
      "ve_pc_up is an AD property (Jacobian flows through density). "
      "ve_dPcup_dsatn is consumed by VEAdvectiveFluxS when capillary=true.");

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
    _dPcup_dsatn(declareProperty<Real>("ve_dPcup_dsatn"))
{
}

void
VEUpscaledCapPressure::computeQpProperties()
{
  const ADReal delta_rho = _density[_qp][1] - _density[_qp][0];
  _pc_up[_qp] = delta_rho * _gravity_magnitude * _h[_qp] + _pc_entry;

  // Sharp-interface formula: d(Pc^up)/d(sat_n) = delta_rho * g * H / (1 - S_wr).
  // Raw value used because density is constant in current fluid property models.
  const Real delta_rho_val = MetaPhysicL::raw_value(delta_rho);
  _dPcup_dsatn[_qp] = delta_rho_val * _gravity_magnitude * _H[_qp] / (1.0 - _S_wr);
}
