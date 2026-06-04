#include "VECapillaryPressureTableUO.h"
#include "DelimitedFileReader.h"
#include "VECsvColumnLookup.h"

#include <algorithm>

registerMooseObject("wombatApp", VECapillaryPressureTableUO);

InputParameters
VECapillaryPressureTableUO::validParams()
{
  InputParameters params = PorousFlowCapillaryPressure::validParams();
  params.addClassDescription(
      "Piecewise-linear capillary pressure UserObject backed by a (Pc, Sw) lookup table. "
      "Intended for use with VEPlumeReconstruction on real field models where the "
      "upscaling workflow (ve_upscale_curves.py) supplies a column-averaged Sw(Pc) table. "
      "Supply the table EITHER inline via pc_points / sw_points OR from a CSV with data_file "
      "(columns Pc, Sw). Use log_extension = false.");

  // Inline table (self-contained input). Optional: omit when using data_file.
  params.addParam<std::vector<Real>>(
      "pc_points",
      "Capillary pressure values [Pa], strictly increasing from 0, defining the "
      "x-axis of the Pc-Sw table. Inline alternative to data_file.");
  params.addParam<std::vector<Real>>(
      "sw_points",
      "True water saturation values [-], strictly decreasing, in one-to-one "
      "correspondence with pc_points.");

  // File table (keeps a large upscaled table out of the input file).
  params.addParam<FileName>(
      "data_file",
      "CSV file with a capillary pressure column and a water-saturation column (e.g. "
      "ve_upscale_curves.py's <prefix>_pc.csv). A header row is auto-detected; when present, "
      "columns are matched by name (pc or Pc; sw or Sw) so their order in the file does not "
      "matter. Without a header the columns are read positionally as Pc, Sw. Alternative to "
      "the inline pc_points / sw_points.");

  // The log extension is designed for parametric curves that diverge at low
  // saturation. A bounded table does not need it and should disable it.
  params.set<bool>("log_extension") = false;

  return params;
}

VECapillaryPressureTableUO::VECapillaryPressureTableUO(const InputParameters & parameters)
  : PorousFlowCapillaryPressure(parameters)
{
  std::vector<Real> pc, sw;

  const bool from_file = isParamValid("data_file");
  const bool from_inline = isParamValid("pc_points") || isParamValid("sw_points");

  if (from_file == from_inline)
    mooseError("VECapillaryPressureTableUO '", name(),
               "': supply EITHER data_file (a CSV with columns Pc, Sw) OR the inline pc_points / "
               "sw_points vectors -- not both, and not neither.");

  if (from_file)
  {
    MooseUtils::DelimitedFileReader reader(getParam<FileName>("data_file"), &_communicator);
    reader.read();
    const auto & cols = reader.getData();
    const auto & names = reader.getNames();

    // If a header is present, match columns by name so the file's column order
    // does not matter. Otherwise (headerless file) fall back to positional order.
    const int pc_col = VECsvColumnLookup::findNamedColumn(names, {"pc", "Pc"});
    const int sw_col = VECsvColumnLookup::findNamedColumn(names, {"sw", "Sw"});
    const bool by_name = (pc_col >= 0) && (sw_col >= 0);

    if (by_name)
    {
      pc = cols[pc_col];
      sw = cols[sw_col];
    }
    else
    {
      if (cols.size() < 2)
        paramError("data_file",
                   "expected at least 2 columns (Pc, Sw) but found ", cols.size(),
                   ". Provide a header row naming the columns (pc or Pc; sw or Sw), or supply "
                   "the two columns in that order without a header.");
      pc = cols[0];
      sw = cols[1];
    }
  }
  else
  {
    if (!isParamValid("pc_points") || !isParamValid("sw_points"))
      mooseError("VECapillaryPressureTableUO '", name(),
                 "': pc_points and sw_points must both be given together (or use data_file).");
    pc = getParam<std::vector<Real>>("pc_points");
    sw = getParam<std::vector<Real>>("sw_points");
  }

  if (pc.size() < 2)
    mooseError("VECapillaryPressureTableUO '", name(), "': at least two table points are required.");
  if (sw.size() != pc.size())
    mooseError("VECapillaryPressureTableUO '", name(),
               "': the Pc and Sw columns must have the same length.");

  for (std::size_t i = 1; i < pc.size(); ++i)
    if (pc[i] <= pc[i - 1])
      mooseError("VECapillaryPressureTableUO '", name(),
                 "': the Pc column must be strictly increasing.");
  for (std::size_t i = 1; i < sw.size(); ++i)
    if (sw[i] >= sw[i - 1])
      mooseError("VECapillaryPressureTableUO '", name(),
                 "': the Sw column must be strictly decreasing.");

  if (pc.front() < 0.0)
    mooseError("VECapillaryPressureTableUO '", name(), "': Pc values must be non-negative.");
  if (sw.back() < _sat_lr)
    mooseError("VECapillaryPressureTableUO '", name(),
               "': the minimum Sw in the table (", sw.back(), ") is below sat_lr (", _sat_lr,
               "). Ensure the table covers down to the residual water saturation.");

  _pc_to_sw.setData(pc, sw);

  // Build the Sw->Pc interpolation with x strictly increasing (reverse both vectors).
  std::vector<Real> sw_rev(sw.rbegin(), sw.rend());
  std::vector<Real> pc_rev(pc.rbegin(), pc.rend());
  _sw_to_pc.setData(sw_rev, pc_rev);
}

Real
VECapillaryPressureTableUO::effectiveSaturation(Real pc, unsigned /*qp*/) const
{
  return effectiveSaturationFromSaturation(_pc_to_sw.sample(pc));
}

Real
VECapillaryPressureTableUO::dEffectiveSaturation(Real pc, unsigned /*qp*/) const
{
  // dSe/dPc = (dSw/dPc) * _dseff_ds   where _dseff_ds = 1 / (1 - sat_lr)
  return _pc_to_sw.sampleDerivative(pc) * _dseff_ds;
}

Real
VECapillaryPressureTableUO::d2EffectiveSaturation(Real /*pc*/, unsigned /*qp*/) const
{
  return 0.0; // piecewise-linear interpolant has zero second derivative
}

Real
VECapillaryPressureTableUO::capillaryPressureCurve(Real saturation, unsigned /*qp*/) const
{
  return _sw_to_pc.sample(saturation);
}

Real
VECapillaryPressureTableUO::dCapillaryPressureCurve(Real saturation, unsigned /*qp*/) const
{
  return _sw_to_pc.sampleDerivative(saturation);
}

Real
VECapillaryPressureTableUO::d2CapillaryPressureCurve(Real /*saturation*/, unsigned /*qp*/) const
{
  return 0.0; // piecewise-linear interpolant has zero second derivative
}
