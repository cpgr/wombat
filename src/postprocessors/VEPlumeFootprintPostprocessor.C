#include "VEPlumeFootprintPostprocessor.h"

registerMooseObject("wombatApp", VEPlumeFootprintPostprocessor);

InputParameters
VEPlumeFootprintPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Integrates the 2-D area where sat_n exceeds a threshold, giving "
      "the CO2 plume footprint (m^2).");

  params.addRequiredCoupledVar("sat_n",
                               "Depth-averaged CO2 saturation (primary variable).");
  params.addParam<Real>("threshold",
                        1.0e-6,
                        "Saturation threshold below which an element is considered "
                        "outside the plume.  Excludes numerical noise at the front.");

  return params;
}

VEPlumeFootprintPostprocessor::VEPlumeFootprintPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _sat_n(coupledValue("sat_n")),
    _threshold(getParam<Real>("threshold"))
{
}

Real
VEPlumeFootprintPostprocessor::computeQpIntegral()
{
  return (_sat_n[_qp] > _threshold) ? 1.0 : 0.0;
}
