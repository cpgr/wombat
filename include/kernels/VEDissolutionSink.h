#pragma once

#include "ADKernelValue.h"

/**
 * VEDissolutionSink
 *
 * Removes free-phase CO2 at the convective dissolution rate ve_dissolution_rate
 * [kg/m^2/s] from VEDissolution. Assign to the CO2 saturation variable (variable=sat_n,
 * the CO2 mass equation). The residual contribution is +ve_dissolution_rate, so that with
 * the storage term d/dt(H phi rho_c S_c) + ... + ve_dissolution_rate = 0 the free CO2 mass
 * decreases (it is a sink). The rate is mass-per-area-per-time, matching the depth-
 * integrated mass storage residual units, so it enters directly with no rho scaling.
 *
 * Off by default in the sense that it is a separate kernel: add it to an input only when
 * dissolution is wanted; existing inputs are unaffected. Jacobian (wrt sat_n, via the
 * gate) is assembled automatically through AD.
 */
class VEDissolutionSink : public ADKernelValue
{
public:
  static InputParameters validParams();
  VEDissolutionSink(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;

  /// Areal CO2 dissolution rate [kg/m^2/s] from VEDissolution.
  const ADMaterialProperty<Real> & _rate;
};
