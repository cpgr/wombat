#include "VERelPermHysteresisUO.h"

#include <algorithm>

registerMooseObject("wombatApp", VERelPermHysteresisUO);

InputParameters
VERelPermHysteresisUO::validParams()
{
  InputParameters params = VERelativePermeability::validParams();
  params.addClassDescription(
      "Hysteretic upscaled relative permeability with residual trapping (saturation-"
      "space sharp-interface Land/Killough). The CO2 curve follows drainage while "
      "sat_n is at its historical max and a linear imbibition scanning curve (to the "
      "Land residual S_gr) once sat_n recedes; the brine curve is non-hysteretic. "
      "Requires the turning point sat_n_max from VESaturationMaxAux (wired through the "
      "VERelPerm / VEFVRelPerm adapters).");

  params.addRangeCheckedParam<Real>(
      "S_wr", 0.0, "S_wr >= 0 & S_wr < 1", "Residual water saturation in the CO2 zone (-).");
  params.addParam<Real>("krn_max", 1.0, "Maximum CO2 (non-wetting) relative permeability (-).");
  params.addParam<Real>("krw_max", 1.0, "Maximum brine (wetting) relative permeability (-).");
  params.addRequiredRangeCheckedParam<Real>(
      "S_gr_max",
      "S_gr_max > 0 & S_gr_max < 1",
      "Maximum residually trapped CO2 saturation after full drainage (-). Must be "
      "< 1 - S_wr (checked in the constructor). Sets the Land coefficient.");

  return params;
}

VERelPermHysteresisUO::VERelPermHysteresisUO(const InputParameters & parameters)
  : VERelativePermeability(parameters),
    _S_wr(getParam<Real>("S_wr")),
    _krn_max(getParam<Real>("krn_max")),
    _krw_max(getParam<Real>("krw_max")),
    _S_gr_max(getParam<Real>("S_gr_max")),
    _C((1.0 - _S_wr) / _S_gr_max - 1.0)
{
  if (_S_gr_max >= 1.0 - _S_wr)
    paramError("S_gr_max",
               "S_gr_max (", _S_gr_max, ") must be < 1 - S_wr (", 1.0 - _S_wr,
               "): the trapped CO2 saturation cannot exceed the movable range.");
}

Real
VERelPermHysteresisUO::effSat(Real sat_n) const
{
  return std::max(0.0, std::min(1.0, sat_n / (1.0 - _S_wr)));
}

Real
VERelPermHysteresisUO::relativePermeability(Real sat_n, unsigned int phase) const
{
  // Two-argument form = primary drainage curve (no history).
  const Real s = effSat(sat_n);
  return (phase == 0) ? _krn_max * s : _krw_max * (1.0 - s);
}

Real
VERelPermHysteresisUO::dRelativePermeability(Real sat_n, unsigned int phase) const
{
  const Real x = sat_n / (1.0 - _S_wr);
  const Real ds = (x > 0.0 && x < 1.0) ? 1.0 / (1.0 - _S_wr) : 0.0;
  return (phase == 0) ? _krn_max * ds : -_krw_max * ds;
}

Real
VERelPermHysteresisUO::relativePermeability(Real sat_n, Real sat_n_max, unsigned int phase) const
{
  // Brine (wetting) is non-hysteretic: depends on total CO2 saturation only.
  if (phase == 1)
    return _krw_max * (1.0 - effSat(sat_n));

  // CO2 (non-wetting): drainage while at/above the historical max, OR when there is
  // no meaningful drainage history (smax ~ 0, e.g. sat_n_max still at its zero IC and
  // a Newton iterate dipping sat_n below it). The denom guard avoids the 0/0 that
  // arises as smax - sgr = C*smax^2/(1+C*smax) -> 0.
  const Real smax = effSat(sat_n_max);
  const Real s = effSat(sat_n);
  const Real sgr = (smax > 0.0) ? smax / (1.0 + _C * smax) : 0.0;
  const Real denom = smax - sgr;
  if (sat_n >= sat_n_max || denom <= 1.0e-12)
    return _krn_max * s; // drainage curve

  const Real mobile = smax * (s - sgr) / denom;
  return _krn_max * std::max(0.0, mobile);
}

Real
VERelPermHysteresisUO::dRelativePermeability(Real sat_n, Real sat_n_max, unsigned int phase) const
{
  if (phase == 1)
  {
    const Real x = sat_n / (1.0 - _S_wr);
    const Real ds = (x > 0.0 && x < 1.0) ? 1.0 / (1.0 - _S_wr) : 0.0;
    return -_krw_max * ds;
  }

  const Real smax = effSat(sat_n_max);
  const Real s = effSat(sat_n);
  const Real sgr = (smax > 0.0) ? smax / (1.0 + _C * smax) : 0.0;
  const Real denom = smax - sgr;
  if (sat_n >= sat_n_max || denom <= 1.0e-12)
    return dRelativePermeability(sat_n, 0); // drainage derivative

  // Nonzero only on the active part of the scanning line, s in (sgr, smax), and
  // inside the saturation clamp. sat_n_max is frozen so smax, sgr carry no derivative.
  if (s > sgr && s < smax)
    return _krn_max * smax / denom / (1.0 - _S_wr);
  return 0.0;
}

Real
VERelPermHysteresisUO::trappedSaturation(Real sat_n, Real sat_n_max) const
{
  const Real smax = effSat(sat_n_max);
  const Real sgr = (smax > 0.0) ? smax / (1.0 + _C * smax) : 0.0;
  const Real denom = smax - sgr;
  if (sat_n >= sat_n_max || denom <= 1.0e-12)
    return 0.0; // drainage / no history: all CO2 connected

  const Real s = effSat(sat_n);
  // Defined so that sat_n - sat_n_trap equals the kr-implied mobile saturation:
  // 0 at the turning point (s = smax), growing to the Land residual (1-S_wr)*sgr at
  // s = sgr. Clamped to [0, sat_n] for safety against round-off / overshoot.
  const Real trapped = (1.0 - _S_wr) * sgr * (smax - s) / denom;
  return std::max(0.0, std::min(trapped, sat_n));
}
