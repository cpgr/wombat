#pragma once

#include "AuxKernel.h"

/**
 * VEPlumeHeightAux
 *
 * Computes the local CO2 plume thickness h [m] from the depth-averaged
 * saturation sat_n and the formation thickness H under the sharp-interface
 * (Nordbotten-Celia) VE assumption:
 *
 *   h = sat_n * H / (1 - S_wr)
 *
 * where S_wr is the residual water saturation in the CO2 zone.  The result
 * is clamped to [0, H] to absorb numerical overshoot at the plume front.
 *
 * Output is an elemental AuxVariable.  Pair with VEPlumeFootprintAux or
 * postprocessors (ElementIntegralVariablePostprocessor, etc.) to obtain
 * aggregate plume statistics.
 *
 * Reads:
 *   - sat_n    : coupled primary variable (depth-averaged CO2 saturation)
 *   - ve_H     : material property from VEGeometry (formation thickness [m])
 *
 * Parameter:
 *   - S_wr     : residual water saturation in the CO2 zone [-] (default 0)
 */
class VEPlumeHeightAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VEPlumeHeightAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  const VariableValue & _sat_n;
  const MaterialProperty<Real> & _H;
  const Real _S_wr;
};
