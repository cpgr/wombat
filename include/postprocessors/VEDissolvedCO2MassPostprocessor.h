#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * VEDissolvedCO2MassPostprocessor
 *
 * Integrates the CO2 mass dissolved in the brine column:
 *
 *   M_diss = integral_Omega  rho_brine * phi_bar * H * (1 - sat_n) * X_co2  dA   [kg]
 *
 * where X_co2 is the dissolved CO2 mass fraction in brine [kg CO2 / kg brine].
 * X_co2 is a coupled AuxVariable written by a future VEDissolution kernel.
 * It defaults to zero when not coupled, making this postprocessor a no-op
 * until dissolution is active.
 *
 * rho_brine is a uniform scalar parameter.
 * Reads ve_H and ve_phi_bar from materials.
 */
class VEDissolvedCO2MassPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();
  VEDissolvedCO2MassPostprocessor(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;

  const VariableValue & _sat_n;
  const VariableValue & _X_co2;
  const MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _phi_bar;
  const Real _rho_brine;
};
