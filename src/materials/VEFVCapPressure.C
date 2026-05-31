#include "VEFVCapPressure.h"

registerMooseObject("wombatApp", VEFVCapPressure);

InputParameters
VEFVCapPressure::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "FV functor ve_pc_up for the upscaled capillary pressure, evaluated at any "
      "face argument so VEFVAdvectiveFlux can compute grad(Pc^up).n via the MOOSE "
      "FV gradient reconstruction. Sharp-interface formula: "
      "Pc^up = delta_rho * |g| * sat_n * (z_top - z_bottom) / (1 - S_wr) + pc_entry. "
      "FE counterpart: VEUpscaledCapPressure.");

  params.addRequiredParam<MooseFunctorName>(
      "sat_n", "Depth-averaged CO2 saturation functor (FV variable name).");
  params.addRequiredParam<MooseFunctorName>(
      "z_top", "Top-surface elevation z_T functor [m].");
  params.addRequiredParam<MooseFunctorName>(
      "z_bottom", "Bottom-surface elevation z_B functor [m].");

  params.addRangeCheckedParam<Real>(
      "S_wr", 0.0, "S_wr >= 0 & S_wr < 1",
      "Residual water saturation [-]. Must match VEPlumeReconstruction and "
      "VEUpscaledCapPressure.");
  params.addRequiredParam<Real>(
      "delta_rho",
      "Density difference rho_brine - rho_co2 [kg/m3]. Positive for buoyant CO2.");
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only the magnitude is used.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Constant capillary entry/fringe pressure offset [Pa].");

  return params;
}

VEFVCapPressure::VEFVCapPressure(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _S_wr(getParam<Real>("S_wr")),
    _delta_rho(getParam<Real>("delta_rho")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _sat_n(getFunctor<ADReal>("sat_n")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom"))
{
  addFunctorProperty<ADReal>(
      "ve_pc_up",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal H = _z_top(r, t) - _z_bottom(r, t);
        const ADReal h = _sat_n(r, t) * H / (1.0 - _S_wr);
        return _delta_rho * _gravity_magnitude * h + _pc_entry;
      });
}
