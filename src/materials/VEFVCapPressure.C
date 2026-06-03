#include "VEFVCapPressure.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("wombatApp", VEFVCapPressure);

InputParameters
VEFVCapPressure::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "FV functors for the upscaled capillary pressure. Sharp-interface formula: "
      "Pc^up = (rho_w - rho_n) * |g| * sat_n * (z_top - z_bottom) / (1 - S_wr) + pc_entry. "
      "The density difference is read per evaluation from the ve_density_n / ve_density_w "
      "functors (VEFVFluidProperties), so it tracks pressure/temperature rather than being "
      "a fixed scalar. Declares ve_pc_up (diagnostic), and the two chain-rule coefficients "
      "ve_dPcup_dsatn = (rho_w-rho_n)*|g|*H/(1-S_wr) and ve_dPcup_dH = "
      "(rho_w-rho_n)*|g|*sat_n/(1-S_wr) that VEFVAdvectiveFlux multiplies by "
      "grad(sat_n).n and grad(H).n to assemble grad(Pc^up).n. "
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
  params.addParam<MooseFunctorName>(
      "density_nw", "ve_density_n",
      "Non-wetting (CO2) density functor [kg/m3], typically from VEFVFluidProperties.");
  params.addParam<MooseFunctorName>(
      "density_w", "ve_density_w",
      "Wetting (brine) density functor [kg/m3], typically from VEFVFluidProperties.");
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
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _sat_n(getFunctor<ADReal>("sat_n")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom")),
    _rho_n(getFunctor<ADReal>("density_nw")),
    _rho_w(getFunctor<ADReal>("density_w"))
{
  // ve_pc_up: the upscaled capillary pressure value (diagnostic / parity with the
  // FE VEUpscaledCapPressure). NOT used for the FV flux gradient -- a material
  // functor's Green-Gauss gradient does not pick up the sat_n Dirichlet BC at a
  // boundary face, so the flux kernel instead multiplies the coefficient below by
  // the sat_n VARIABLE gradient (gradUDotNormal), which is boundary-aware.
  // delta_rho = rho_w - rho_n is read from the density functors at each evaluation,
  // so it follows the EOS (pressure/temperature) rather than being a fixed scalar.
  addFunctorProperty<ADReal>(
      "ve_pc_up",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
        const ADReal H = _z_top(r, t) - _z_bottom(r, t);
        const ADReal h = _sat_n(r, t) * H / (1.0 - _S_wr);
        return delta_rho * _gravity_magnitude * h + _pc_entry;
      });

  // ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = (rho_w - rho_n) * |g| * H / (1 - S_wr).
  // Consumed by VEFVAdvectiveFlux (capillary=true): grad(Pc^up).n is formed as
  // ve_dPcup_dsatn(face) * grad(sat_n).n. Only its VALUE is evaluated (no gradient
  // of a material functor is taken), so the BC-awareness comes from grad(sat_n).
  // The density's pressure derivative is stripped here (raw_value) -- as in the FE
  // VEUpscaledCapPressure coefficient -- because the sat_n Jacobian flows through
  // grad(sat_n) and the density-vs-pp_top sensitivity of this coefficient is a
  // second-order effect (and exactly zero for constant fluid properties).
  addFunctorProperty<ADReal>(
      "ve_dPcup_dsatn",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const Real delta_rho = MetaPhysicL::raw_value(_rho_w(r, t) - _rho_n(r, t));
        const Real H = MetaPhysicL::raw_value(_z_top(r, t) - _z_bottom(r, t));
        return delta_rho * _gravity_magnitude * H / (1.0 - _S_wr);
      });

  // ve_dPcup_dH = d(Pc^up)/d(H) = (rho_w - rho_n) * |g| * sat_n / (1 - S_wr). The
  // second half of the chain rule grad(Pc^up).n = ve_dPcup_dsatn*grad(sat_n).n
  // + ve_dPcup_dH*grad(H).n for laterally varying thickness. VEFVAdvectiveFlux
  // multiplies its face VALUE by grad(H).n formed from the z_top/z_bottom functor
  // gradients (fixed geometry, BC-safe), so a Dirichlet sat_n inlet still enters
  // through the boundary-aware sat_n VALUE here. Kept AD (as in FE ve_dPcup_dH) so
  // the d/d(sat_n) derivative rides on this coefficient.
  addFunctorProperty<ADReal>(
      "ve_dPcup_dH",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
        return delta_rho * _gravity_magnitude * _sat_n(r, t) / (1.0 - _S_wr);
      });
}
