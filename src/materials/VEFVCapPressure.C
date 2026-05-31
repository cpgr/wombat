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
  // ve_pc_up: the upscaled capillary pressure value (diagnostic / parity with the
  // FE VEUpscaledCapPressure). NOT used for the FV flux gradient -- a material
  // functor's Green-Gauss gradient does not pick up the sat_n Dirichlet BC at a
  // boundary face, so the flux kernel instead multiplies the coefficient below by
  // the sat_n VARIABLE gradient (gradUDotNormal), which is boundary-aware.
  addFunctorProperty<ADReal>(
      "ve_pc_up",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal H = _z_top(r, t) - _z_bottom(r, t);
        const ADReal h = _sat_n(r, t) * H / (1.0 - _S_wr);
        return _delta_rho * _gravity_magnitude * h + _pc_entry;
      });

  // ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = delta_rho * |g| * H / (1 - S_wr).
  // Consumed by VEFVAdvectiveFlux (capillary=true): grad(Pc^up).n is formed as
  // ve_dPcup_dsatn(face) * grad(sat_n).n. Only its VALUE is evaluated (no gradient
  // of a material functor is taken), so the BC-awareness comes from grad(sat_n).
  // Mirrors the FE VEUpscaledCapPressure ve_dPcup_dsatn coefficient.
  addFunctorProperty<ADReal>(
      "ve_dPcup_dsatn",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal H = _z_top(r, t) - _z_bottom(r, t);
        return _delta_rho * _gravity_magnitude * H / (1.0 - _S_wr);
      });
}
