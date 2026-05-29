#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * VEPlumeFootprintPostprocessor
 *
 * Integrates the area of elements where the depth-averaged CO2 saturation
 * sat_n exceeds a small threshold, giving the 2-D footprint of the plume:
 *
 *   A_plume = integral_Omega  I(sat_n > threshold)  dA   [m^2]
 *
 * The threshold defaults to 1e-6 to exclude numerical noise at the front.
 * Pair with VEMobileCO2MassPostprocessor and VETrappedCO2MassPostprocessor
 * for a full mass-balance summary.
 */
class VEPlumeFootprintPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();
  VEPlumeFootprintPostprocessor(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;

  const VariableValue & _sat_n;
  const Real _threshold;
};
