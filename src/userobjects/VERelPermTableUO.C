#include "VERelPermTableUO.h"

registerMooseObject("wombatApp", VERelPermTableUO);

InputParameters
VERelPermTableUO::validParams()
{
  InputParameters params = VERelativePermeability::validParams();
  params.addClassDescription(
      "Tabulated upscaled relative permeability: piecewise-linear lookup of "
      "kr_n(sat_n) and kr_w(sat_n) from (sat_n, kr) tables, e.g. produced by the "
      "fine-scale upscaling workflow. Field-case counterpart of VERelPermSharpUO.");

  params.addRequiredParam<std::vector<Real>>(
      "sat_n_points",
      "Depth-averaged CO2 saturation sample points [-], strictly increasing, "
      "defining the x-axis shared by both relperm tables.");
  params.addRequiredParam<std::vector<Real>>(
      "krn_points", "CO2 (non-wetting) relative permeability at sat_n_points [-].");
  params.addRequiredParam<std::vector<Real>>(
      "krw_points", "Brine (wetting) relative permeability at sat_n_points [-].");

  return params;
}

VERelPermTableUO::VERelPermTableUO(const InputParameters & parameters)
  : VERelativePermeability(parameters)
{
  const auto & sat = getParam<std::vector<Real>>("sat_n_points");
  const auto & krn = getParam<std::vector<Real>>("krn_points");
  const auto & krw = getParam<std::vector<Real>>("krw_points");

  if (sat.size() < 2)
    paramError("sat_n_points", "At least two points are required.");
  if (krn.size() != sat.size())
    paramError("krn_points", "krn_points and sat_n_points must have the same length.");
  if (krw.size() != sat.size())
    paramError("krw_points", "krw_points and sat_n_points must have the same length.");

  for (std::size_t i = 1; i < sat.size(); ++i)
    if (sat[i] <= sat[i - 1])
      paramError("sat_n_points", "Values must be strictly increasing.");

  _krn.setData(sat, krn);
  _krw.setData(sat, krw);
}

Real
VERelPermTableUO::relativePermeability(Real sat_n, unsigned int phase) const
{
  return (phase == 0) ? _krn.sample(sat_n) : _krw.sample(sat_n);
}

Real
VERelPermTableUO::dRelativePermeability(Real sat_n, unsigned int phase) const
{
  return (phase == 0) ? _krn.sampleDerivative(sat_n) : _krw.sampleDerivative(sat_n);
}
