#include "VERelPermSharp.h"

registerMooseObject("wombatApp", VERelPermSharp);

InputParameters
VERelPermSharp::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Sharp-interface upscaled relative permeabilities (Nordbotten-Celia). "
      "kr_n_up = krn_max * sat_n / (1 - S_wr), "
      "kr_w_up = krw_max * (1 - sat_n / (1 - S_wr)). "
      "Reads ve_saturation, writes ve_relperm.");

  params.addParam<Real>("S_wr",    0.0, "Residual water saturation in the CO2 zone [-].");
  params.addParam<Real>("krn_max", 1.0, "Maximum CO2 (non-wetting) relative permeability [-].");
  params.addParam<Real>("krw_max", 1.0, "Maximum brine (wetting) relative permeability [-].");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VERelPermSharp::VERelPermSharp(const InputParameters & parameters)
  : Material(parameters),
    _S_wr(getParam<Real>("S_wr")),
    _krn_max(getParam<Real>("krn_max")),
    _krw_max(getParam<Real>("krw_max")),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _relperm(declareADProperty<std::vector<Real>>("ve_relperm"))
{
  if (_S_wr < 0.0 || _S_wr >= 1.0)
    paramError("S_wr", "Residual water saturation must be in [0, 1).");
}

void
VERelPermSharp::computeQpProperties()
{
  _relperm[_qp].resize(2);

  // Normalised CO2 saturation, clamped to [0,1] to absorb numerical overshoot.
  // AD propagates derivatives of sat_n through the clamp automatically.
  const ADReal eff_sat =
      std::max(ADReal(0.0), std::min(ADReal(1.0), _saturation[_qp][0] / (1.0 - _S_wr)));

  _relperm[_qp][0] = _krn_max * eff_sat;          // kr_n_up
  _relperm[_qp][1] = _krw_max * (1.0 - eff_sat);  // kr_w_up
}
