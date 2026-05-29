#include "VECapillaryPressureTable.h"

#include <algorithm>

registerMooseObject("wombatApp", VECapillaryPressureTable);

InputParameters
VECapillaryPressureTable::validParams()
{
  InputParameters params = PorousFlowCapillaryPressure::validParams();
  params.addClassDescription(
      "Piecewise-linear capillary pressure UserObject backed by a (Pc, Sw) lookup table. "
      "Intended for use with VEPlumeReconstruction on real field models where the "
      "upscaling workflow supplies a column-averaged Sw(Pc) table that encodes all "
      "vertical facies heterogeneity. Use log_extension = false.");

  params.addRequiredParam<std::vector<Real>>(
      "pc_points",
      "Capillary pressure values [Pa], strictly increasing from 0, defining the "
      "x-axis of the Pc-Sw table (output from the upscaling workflow).");
  params.addRequiredParam<std::vector<Real>>(
      "sw_points",
      "True water saturation values [-], strictly decreasing, in one-to-one "
      "correspondence with pc_points (output from the upscaling workflow).");

  // The log extension is designed for parametric curves that diverge at low
  // saturation. A bounded table does not need it and should disable it.
  params.set<bool>("log_extension") = false;

  return params;
}

VECapillaryPressureTable::VECapillaryPressureTable(const InputParameters & parameters)
  : PorousFlowCapillaryPressure(parameters)
{
  const auto & pc = getParam<std::vector<Real>>("pc_points");
  const auto & sw = getParam<std::vector<Real>>("sw_points");

  if (pc.size() < 2)
    paramError("pc_points", "At least two points are required.");
  if (sw.size() != pc.size())
    paramError("sw_points", "sw_points and pc_points must have the same length.");

  for (std::size_t i = 1; i < pc.size(); ++i)
    if (pc[i] <= pc[i - 1])
      paramError("pc_points", "Values must be strictly increasing.");
  for (std::size_t i = 1; i < sw.size(); ++i)
    if (sw[i] >= sw[i - 1])
      paramError("sw_points", "Values must be strictly decreasing.");

  if (pc.front() < 0.0)
    paramError("pc_points", "Capillary pressure values must be non-negative.");
  if (sw.back() < _sat_lr)
    paramError("sw_points",
               "Minimum Sw in the table is below sat_lr. "
               "Ensure the table covers down to the residual water saturation.");

  _pc_to_sw.setData(pc, sw);

  // Build the Sw→Pc interpolation with x strictly increasing (reverse both vectors).
  std::vector<Real> sw_rev(sw.rbegin(), sw.rend());
  std::vector<Real> pc_rev(pc.rbegin(), pc.rend());
  _sw_to_pc.setData(sw_rev, pc_rev);
}

Real
VECapillaryPressureTable::effectiveSaturation(Real pc, unsigned /*qp*/) const
{
  return effectiveSaturationFromSaturation(_pc_to_sw.sample(pc));
}

Real
VECapillaryPressureTable::dEffectiveSaturation(Real pc, unsigned /*qp*/) const
{
  // dSe/dPc = (dSw/dPc) * _dseff_ds   where _dseff_ds = 1 / (1 - sat_lr)
  return _pc_to_sw.sampleDerivative(pc) * _dseff_ds;
}

Real
VECapillaryPressureTable::d2EffectiveSaturation(Real /*pc*/, unsigned /*qp*/) const
{
  return 0.0; // piecewise-linear interpolant has zero second derivative
}

Real
VECapillaryPressureTable::capillaryPressureCurve(Real saturation, unsigned /*qp*/) const
{
  return _sw_to_pc.sample(saturation);
}

Real
VECapillaryPressureTable::dCapillaryPressureCurve(Real saturation, unsigned /*qp*/) const
{
  return _sw_to_pc.sampleDerivative(saturation);
}

Real
VECapillaryPressureTable::d2CapillaryPressureCurve(Real /*saturation*/, unsigned /*qp*/) const
{
  return 0.0; // piecewise-linear interpolant has zero second derivative
}
