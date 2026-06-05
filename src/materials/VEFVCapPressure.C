#include "VEFVCapPressure.h"
#include "PorousFlowCapillaryPressure.h"

#include "metaphysicl/raw_type.h"

#include <algorithm>

registerMooseObject("wombatApp", VEFVCapPressure);

InputParameters
VEFVCapPressure::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "FV functors for the upscaled capillary pressure. grad(Pc^up) = delta_rho*|g|*grad(h) "
      "in both modes, so the coefficients are the chain-rule partials of h(sat_n, H): "
      "sharp_interface uses h = sat_n*(z_top-z_bottom)/(1-S_wr); capillary_fringe Newton-inverts "
      "the column equilibrium with a PorousFlowCapillaryPressure UO and uses S_n(h) = "
      "1 - Sw(delta_rho*|g|*h). The density difference is read per evaluation from the "
      "ve_density_n / ve_density_w functors. Declares ve_pc_up (diagnostic) and the two "
      "chain-rule coefficients ve_dPcup_dsatn, ve_dPcup_dH that VEFVAdvectiveFlux multiplies "
      "by grad(sat_n).n and grad(H).n. FE counterpart: VEUpscaledCapPressure.");

  params.addRequiredParam<MooseFunctorName>(
      "sat_n", "Depth-averaged CO2 saturation functor (FV variable name).");
  params.addRequiredParam<MooseFunctorName>(
      "z_top", "Top-surface elevation z_T functor (m).");
  params.addRequiredParam<MooseFunctorName>(
      "z_bottom", "Bottom-surface elevation z_B functor (m).");

  MooseEnum modes("sharp_interface capillary_fringe", "sharp_interface");
  params.addParam<MooseEnum>(
      "mode",
      modes,
      "Reconstruction mode. Must match VEPlumeReconstruction's mode where one is used.");
  params.addParam<UserObjectName>(
      "pc_uo",
      "PorousFlowCapillaryPressure UserObject defining the Sw(Pc) curve. Required in "
      "capillary_fringe mode.");

  params.addRangeCheckedParam<Real>(
      "S_wr", 0.0, "S_wr >= 0 & S_wr < 1",
      "Residual water saturation (-). Must match VEPlumeReconstruction and "
      "VEUpscaledCapPressure.");
  params.addParam<MooseFunctorName>(
      "density_nw", "ve_density_n",
      "Non-wetting (CO2) density functor (kg/m3), typically from VEFVFluidProperties.");
  params.addParam<MooseFunctorName>(
      "density_w", "ve_density_w",
      "Wetting (brine) density functor (kg/m3), typically from VEFVFluidProperties.");
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector (m/s2). Only the magnitude is used.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Constant capillary entry/fringe pressure offset (Pa).");

  return params;
}

VEFVCapPressure::VEFVCapPressure(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _mode(getParam<MooseEnum>("mode")),
    _S_wr(getParam<Real>("S_wr")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _pc_entry(getParam<Real>("pc_entry")),
    _sat_n(getFunctor<ADReal>("sat_n")),
    _z_top(getFunctor<ADReal>("z_top")),
    _z_bottom(getFunctor<ADReal>("z_bottom")),
    _rho_n(getFunctor<ADReal>("density_nw")),
    _rho_w(getFunctor<ADReal>("density_w")),
    _pc_uo(nullptr)
{
  if (_mode == "capillary_fringe")
  {
    if (!isParamValid("pc_uo"))
      paramError("pc_uo", "Must be supplied for capillary_fringe mode.");
    _pc_uo = &getUserObject<PorousFlowCapillaryPressure>("pc_uo");
  }

  // ve_pc_up: the upscaled capillary pressure VALUE (diagnostic / parity with the
  // FE VEUpscaledCapPressure). NOT used for the FV flux gradient -- a material
  // functor's Green-Gauss gradient does not pick up the sat_n Dirichlet BC at a
  // boundary face, so the flux kernel multiplies the coefficients below by the
  // sat_n VARIABLE gradient (gradUDotNormal), which is boundary-aware.
  addFunctorProperty<ADReal>(
      "ve_pc_up",
      [this](const auto & r, const auto & t) -> ADReal
      {
        const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
        const ADReal H = _z_top(r, t) - _z_bottom(r, t);
        const ADReal h = plumeThickness(delta_rho, H, _sat_n(r, t));
        return delta_rho * _gravity_magnitude * h + _pc_entry;
      });

  // ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = delta_rho*|g|*dh/dsat_n. Consumed by
  // VEFVAdvectiveFlux (capillary != none): grad(Pc^up).n = ve_dPcup_dsatn(face) *
  // grad(sat_n).n. The functor VALUE carries d/d(sat_n) (sharp: zero; fringe: via
  // the AD seed on h and S_n), and it multiplies the AD gradUDotNormal, so the
  // flux Jacobian is AD-consistent in both modes.
  if (_mode == "sharp_interface")
    addFunctorProperty<ADReal>(
        "ve_dPcup_dsatn",
        [this](const auto & r, const auto & t) -> ADReal
        {
          // sharp dh/dsat_n = H/(1-S_wr), constant in sat_n. delta_rho stripped
          // to raw (as in FE sharp): the sat_n Jacobian rides on grad(sat_n) and
          // the density-vs-pressure sensitivity is second order (zero for
          // constant fluid properties).
          const Real delta_rho = MetaPhysicL::raw_value(_rho_w(r, t) - _rho_n(r, t));
          const Real H = MetaPhysicL::raw_value(_z_top(r, t) - _z_bottom(r, t));
          return delta_rho * _gravity_magnitude * H / (1.0 - _S_wr);
        });
  else
    addFunctorProperty<ADReal>(
        "ve_dPcup_dsatn",
        [this](const auto & r, const auto & t) -> ADReal
        {
          const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
          const ADReal H = _z_top(r, t) - _z_bottom(r, t);
          const ADReal h = plumeThickness(delta_rho, H, _sat_n(r, t));
          return delta_rho * _gravity_magnitude * H / saturationTop(delta_rho, h);
        });

  // ve_dPcup_dH = d(Pc^up)/d(H) = delta_rho*|g|*dh/dH. The second half of the
  // chain rule grad(Pc^up).n = ve_dPcup_dsatn*grad(sat_n).n + ve_dPcup_dH*grad(H).n
  // for laterally varying thickness. Kept AD so the d/d(sat_n) derivative rides on
  // this coefficient (grad(H).n is fixed geometry).
  if (_mode == "sharp_interface")
    addFunctorProperty<ADReal>(
        "ve_dPcup_dH",
        [this](const auto & r, const auto & t) -> ADReal
        {
          // sharp dh/dH = sat_n/(1-S_wr).
          const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
          return delta_rho * _gravity_magnitude * _sat_n(r, t) / (1.0 - _S_wr);
        });
  else
    addFunctorProperty<ADReal>(
        "ve_dPcup_dH",
        [this](const auto & r, const auto & t) -> ADReal
        {
          // fringe dh/dH = sat_n/S_n(h).
          const ADReal delta_rho = _rho_w(r, t) - _rho_n(r, t);
          const ADReal H = _z_top(r, t) - _z_bottom(r, t);
          const ADReal sat_n = _sat_n(r, t);
          const ADReal h = plumeThickness(delta_rho, H, sat_n);
          return delta_rho * _gravity_magnitude * sat_n / saturationTop(delta_rho, h);
        });
}

Real
VEFVCapPressure::satNBar(Real h, Real H, Real delta_rho) const
{
  // Trapezoidal integration of S_n(zeta) = 1 - Sw(Pc(zeta)) over [0, h], with
  // Pc(zeta) = delta_rho*|g|*zeta. Mirrors VEPlumeReconstruction::satNBar.
  const unsigned int N = 64;
  const Real dz = h / N;
  Real sum = 0.0;
  for (unsigned int i = 0; i <= N; ++i)
  {
    const Real pc = delta_rho * _gravity_magnitude * i * dz;
    const Real Sn = 1.0 - _pc_uo->saturation(pc);
    sum += (i == 0 || i == N) ? 0.5 * Sn : Sn;
  }
  return sum * dz / H;
}

ADReal
VEFVCapPressure::plumeThickness(const ADReal & delta_rho,
                                const ADReal & H,
                                const ADReal & sat_n) const
{
  if (_mode == "sharp_interface")
    return sat_n * H / (1.0 - _S_wr);

  // capillary_fringe: Newton inversion on raw values (mirrors VEPlumeReconstruction).
  const Real drho = MetaPhysicL::raw_value(delta_rho);
  const Real H_val = MetaPhysicL::raw_value(H);
  const Real raw_sat = MetaPhysicL::raw_value(sat_n);

  Real h_val = std::max(0.0, std::min(raw_sat * H_val / (1.0 - _S_wr), H_val)); // warm start

  for (unsigned int iter = 0; iter < 20; ++iter)
  {
    const Real F = satNBar(h_val, H_val, drho) - raw_sat;
    if (std::abs(F) < 1.0e-12)
      break;
    const Real Sn_top = 1.0 - _pc_uo->saturation(drho * _gravity_magnitude * h_val);
    if (Sn_top < 1.0e-14)
      break;
    h_val -= F * H_val / Sn_top;
    h_val = std::max(0.0, std::min(h_val, H_val));
  }

  if (h_val <= 0.0 || h_val >= H_val)
    return ADReal(h_val);

  // Linearise: h(sat_n) ~= h_val + (H/S_n(h)) * (sat_n - raw_sat).
  const Real Sn_h = 1.0 - _pc_uo->saturation(drho * _gravity_magnitude * h_val);
  const Real dh_dsatn = (Sn_h > 1.0e-14) ? H_val / Sn_h : H_val / 1.0e-14;
  return h_val + dh_dsatn * (sat_n - raw_sat);
}

ADReal
VEFVCapPressure::saturationTop(const ADReal & delta_rho, const ADReal & h) const
{
  // S_n(h) = 1 - Sw(delta_rho*|g|*h). The AD saturation overload propagates h's
  // d/d(sat_n) through dSw/dPc. Floor near zero (plume nose) to keep H/S_n finite.
  const ADReal Sn = 1.0 - _pc_uo->saturation(delta_rho * _gravity_magnitude * h);
  if (MetaPhysicL::raw_value(Sn) < 1.0e-12)
    return ADReal(1.0e-12);
  return Sn;
}
