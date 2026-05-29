#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * VETrappedCO2MassPostprocessor
 *
 * Integrates the depth-averaged residually trapped CO2 mass:
 *
 *   M_trap = integral_Omega  rho_co2 * phi_bar * H * sat_n_trap  dA   [kg]
 *
 * sat_n_trap is a coupled AuxVariable written by the VEHysteresis material
 * (stateful, tracks the scanning-curve trapped saturation).  It defaults to
 * zero when not coupled, making this postprocessor a no-op until hysteresis
 * is active -- which is the correct behaviour for v1 sharp-interface cases.
 *
 * rho_co2 is a uniform scalar parameter.
 * Reads ve_H and ve_phi_bar from materials.
 */
class VETrappedCO2MassPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();
  VETrappedCO2MassPostprocessor(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;

  const VariableValue & _sat_n_trap;
  const MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _phi_bar;
  const Real _rho_co2;
};
