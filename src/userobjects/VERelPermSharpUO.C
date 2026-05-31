#include "VERelPermSharpUO.h"

#include <algorithm>

registerMooseObject("wombatApp", VERelPermSharpUO);

InputParameters
VERelPermSharpUO::validParams()
{
  InputParameters params = VERelativePermeability::validParams();
  params.addClassDescription(
      "Sharp-interface (Nordbotten-Celia) upscaled relative permeability: "
      "kr_n = krn_max * s_eff and kr_w = krw_max * (1 - s_eff) with "
      "s_eff = clamp(sat_n / (1 - S_wr), 0, 1).");
  params.addRangeCheckedParam<Real>("S_wr",
                                    0.0,
                                    "S_wr >= 0 & S_wr < 1",
                                    "Residual water saturation in the CO2 zone [-].");
  params.addParam<Real>("krn_max", 1.0, "Maximum CO2 (non-wetting) relative permeability [-].");
  params.addParam<Real>("krw_max", 1.0, "Maximum brine (wetting) relative permeability [-].");
  return params;
}

VERelPermSharpUO::VERelPermSharpUO(const InputParameters & parameters)
  : VERelativePermeability(parameters),
    _S_wr(getParam<Real>("S_wr")),
    _krn_max(getParam<Real>("krn_max")),
    _krw_max(getParam<Real>("krw_max"))
{
}

Real
VERelPermSharpUO::relativePermeability(Real sat_n, unsigned int phase) const
{
  const Real eff = std::max(0.0, std::min(1.0, sat_n / (1.0 - _S_wr)));
  return (phase == 0) ? _krn_max * eff : _krw_max * (1.0 - eff);
}

Real
VERelPermSharpUO::dRelativePermeability(Real sat_n, unsigned int phase) const
{
  // d(s_eff)/d(sat_n) is 1/(1 - S_wr) on the open interval and 0 at/outside the
  // clamp. The strict inequalities reproduce the AD-through-clamp behaviour of
  // the previous material exactly (std::max/std::min return the constant first
  // argument on a tie, i.e. zero derivative at the kinks).
  const Real x = sat_n / (1.0 - _S_wr);
  const Real deff = (x > 0.0 && x < 1.0) ? 1.0 / (1.0 - _S_wr) : 0.0;
  return (phase == 0) ? _krn_max * deff : -_krw_max * deff;
}
