#pragma once

#include "FVElementalKernel.h"

/**
 * VEFVDissolutionSink
 *
 * Finite-volume analogue of VEDissolutionSink. Removes free-phase CO2 at the convective
 * dissolution rate ve_dissolution_rate [kg/m^2/s] (from VEDissolution) on the CO2 mass
 * equation (variable=sat_n). The residual contribution is +ve_dissolution_rate, the same
 * sign convention as the FV mass storage term VEFVMassTimeDerivative (elemental terms share
 * the FE sign; only the flux kernels differ in the div(F) sign), so the free CO2 mass
 * decreases. The rate is an areal mass flux, matching the depth-integrated storage units,
 * so it enters directly with no rho scaling. Jacobian (wrt sat_n via the gate) via AD.
 */
class VEFVDissolutionSink : public FVElementalKernel
{
public:
  static InputParameters validParams();
  VEFVDissolutionSink(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  /// Areal CO2 dissolution rate [kg/m^2/s] from VEDissolution.
  const ADMaterialProperty<Real> & _rate;
};
