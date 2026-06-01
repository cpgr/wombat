#include "VEDissolution.h"

registerMooseObject("wombatApp", VEDissolution);

InputParameters
VEDissolution::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Convective CO2 dissolution rate ve_dissolution_rate = q0 * gate(sat_n) * "
      "capacity(c_diss) [kg/m^2/s], an areal mass flux consumed by VEDissolutionSink and "
      "accumulated by VEDissolvedCO2Aux. q0 is the constant-flux convective rate; the gate "
      "rate-limits as the plume thins; the optional capacity factor stops dissolution at "
      "column CO2 saturation. Sink-only (Option A): dissolved CO2 is an immobile inventory.");

  params.addRequiredRangeCheckedParam<Real>(
      "dissolution_flux",
      "dissolution_flux >= 0",
      "Constant-flux convective dissolution rate q0 per unit map-view area [kg/m^2/s].");
  params.addRangeCheckedParam<Real>(
      "s_ref",
      0.05,
      "s_ref > 0",
      "Gate reference CO2 saturation [-]: full dissolution for sat_n >= s_ref, tapering "
      "linearly to zero as sat_n -> 0 so the sink cannot over-deplete the free phase.");
  params.addCoupledVar(
      "dissolved_co2",
      "Optional lagged areal dissolved CO2 mass c_diss [kg/m^2] (from VEDissolvedCO2Aux). "
      "Supply together with c_cap to cap dissolution at the column CO2 capacity.");
  params.addRangeCheckedParam<Real>(
      "c_cap",
      "c_cap > 0",
      "Column CO2 capacity [kg/m^2] (e.g. rho_brine*phi*H*X_sat). Dissolution stops as "
      "c_diss -> c_cap. Requires dissolved_co2; omit both to disable the cap.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEDissolution::VEDissolution(const InputParameters & parameters)
  : Material(parameters),
    _q0(getParam<Real>("dissolution_flux")),
    _s_ref(getParam<Real>("s_ref")),
    _has_cap(isCoupled("dissolved_co2") && isParamValid("c_cap")),
    _c_cap(_has_cap ? getParam<Real>("c_cap") : 1.0),
    _c_diss(isCoupled("dissolved_co2") ? coupledValue("dissolved_co2") : _zero),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _rate(declareADProperty<Real>("ve_dissolution_rate"))
{
  if (isCoupled("dissolved_co2") != isParamValid("c_cap"))
    paramError("c_cap",
               "Provide both dissolved_co2 and c_cap to enable the capacity cap, or "
               "neither to disable it.");
}

void
VEDissolution::computeQpProperties()
{
  const ADReal sat_n = _saturation[_qp][0];

  // Plume-presence gate = min(1, sat_n / s_ref): exactly 1 (and AD-flat) where the plume
  // is present, linear taper to 0 as sat_n -> 0 so the sink cannot drive sat_n negative.
  ADReal gate;
  if (MetaPhysicL::raw_value(sat_n) >= _s_ref)
    gate = 1.0;
  else
    gate = sat_n / _s_ref;

  // Optional column-capacity factor (lagged c_diss -> frozen, not AD).
  Real cap = 1.0;
  if (_has_cap)
    cap = std::max(0.0, 1.0 - _c_diss[_qp] / _c_cap);

  _rate[_qp] = _q0 * gate * cap;
}
