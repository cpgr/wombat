#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * VEMobileCO2MassPostprocessor
 *
 * Integrates the depth-averaged mobile CO2 mass over the 2-D domain:
 *
 *   M_mobile = integral_Omega  rho_co2 * phi_bar * H * sat_n_mobile  dA   [kg]
 *
 * where sat_n_mobile = sat_n - sat_n_trap (free CO2 saturation).
 * Trapped saturation sat_n_trap may optionally be coupled as an AuxVariable
 * (produced by a future VEHysteresis material).  When not coupled it defaults
 * to zero, so the full sat_n is treated as mobile -- correct for the
 * no-hysteresis verification cases.
 *
 * rho_co2 is a uniform scalar parameter (valid where the CO2 density is constant,
 * i.e. VEFluidProperties driven by ConstantFluidProperties).
 * Reads ve_H and ve_phi_bar from materials.
 */
class VEMobileCO2MassPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();
  VEMobileCO2MassPostprocessor(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;

  const VariableValue & _sat_n;
  const VariableValue & _sat_n_trap;
  const MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _phi_bar;
  const Real _rho_co2;
};
