#include "VERelPermTableUO.h"
#include "DelimitedFileReader.h"

registerMooseObject("wombatApp", VERelPermTableUO);

InputParameters
VERelPermTableUO::validParams()
{
  InputParameters params = VERelativePermeability::validParams();
  params.addClassDescription(
      "Tabulated upscaled relative permeability: piecewise-linear lookup of "
      "kr_n(sat_n) and kr_w(sat_n) from (sat_n, kr) tables, e.g. produced by the "
      "fine-scale upscaling workflow (ve_upscale_curves.py). Field-case counterpart of "
      "VERelPermSharpUO. Supply the table EITHER inline via sat_n_points / krn_points / "
      "krw_points OR from a CSV with data_file (columns sat_n, kr_n, kr_w).");

  // Inline table (self-contained input). Optional: omit when using data_file.
  params.addParam<std::vector<Real>>(
      "sat_n_points",
      "Depth-averaged CO2 saturation sample points [-], strictly increasing, "
      "defining the x-axis shared by both relperm tables. Inline alternative to data_file.");
  params.addParam<std::vector<Real>>(
      "krn_points", "CO2 (non-wetting) relative permeability at sat_n_points [-].");
  params.addParam<std::vector<Real>>(
      "krw_points", "Brine (wetting) relative permeability at sat_n_points [-].");

  // File table (keeps a large upscaled table out of the input file).
  params.addParam<FileName>(
      "data_file",
      "CSV file with three columns sat_n, kr_n, kr_w (e.g. ve_upscale_curves.py's "
      "<prefix>_relperm.csv). A header row is auto-detected. Alternative to the inline "
      "sat_n_points / krn_points / krw_points.");

  return params;
}

VERelPermTableUO::VERelPermTableUO(const InputParameters & parameters)
  : VERelativePermeability(parameters)
{
  std::vector<Real> sat, krn, krw;

  const bool from_file = isParamValid("data_file");
  const bool from_inline =
      isParamValid("sat_n_points") || isParamValid("krn_points") || isParamValid("krw_points");

  if (from_file == from_inline)
    mooseError("VERelPermTableUO '", name(),
               "': supply EITHER data_file (a CSV with columns sat_n, kr_n, kr_w) OR the inline "
               "sat_n_points / krn_points / krw_points vectors -- not both, and not neither.");

  if (from_file)
  {
    MooseUtils::DelimitedFileReader reader(getParam<FileName>("data_file"), &_communicator);
    reader.read();
    const auto & cols = reader.getData();
    if (cols.size() < 3)
      paramError("data_file",
                 "expected at least 3 columns (sat_n, kr_n, kr_w) but found ", cols.size(), ".");
    sat = cols[0];
    krn = cols[1];
    krw = cols[2];
  }
  else
  {
    if (!isParamValid("sat_n_points") || !isParamValid("krn_points") ||
        !isParamValid("krw_points"))
      mooseError("VERelPermTableUO '", name(),
                 "': sat_n_points, krn_points and krw_points must all be given together (or use "
                 "data_file instead).");
    sat = getParam<std::vector<Real>>("sat_n_points");
    krn = getParam<std::vector<Real>>("krn_points");
    krw = getParam<std::vector<Real>>("krw_points");
  }

  if (sat.size() < 2)
    mooseError("VERelPermTableUO '", name(), "': at least two table points are required.");
  if (krn.size() != sat.size() || krw.size() != sat.size())
    mooseError("VERelPermTableUO '", name(),
               "': the sat_n, kr_n and kr_w columns must have the same length.");

  for (std::size_t i = 1; i < sat.size(); ++i)
    if (sat[i] <= sat[i - 1])
      mooseError("VERelPermTableUO '", name(),
                 "': the sat_n column must be strictly increasing.");

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
