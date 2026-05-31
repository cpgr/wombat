#pragma once

#include "VEAdvectiveFlux.h"

/**
 * VEAdvectiveFluxS
 *
 * Depth-integrated Darcy advective flux kernel for use on the saturation
 * variable (variable = sat_n).
 *
 * Because the kernel variable is sat_n, _grad_u[_qp] holds
 * grad(sat_n), which is NOT the pressure gradient needed for Darcy
 * flow. Instead, this kernel couples explicitly to pp_top and uses
 * adCoupledGradient("pp_top") to obtain grad(pp_top) with full AD
 * derivatives (so off-diagonal Jacobian entries wrt pp_top are correct).
 *
 * This variant is used for the CO2 (non-wetting phase) mass conservation
 * equation:
 *
 *   fluid_phase = 0  (CO2),  variable = sat_n,  coupled pp_top
 */
class VEAdvectiveFluxS : public VEAdvectiveFluxBase
{
public:
  static InputParameters validParams();
  VEAdvectiveFluxS(const InputParameters & parameters);

protected:
  ADRealVectorValue precomputeQpResidual() override;

  /// AD gradient of pp_top -- carries off-diagonal Jacobian entries wrt pressure.
  const ADVariableGradient & _grad_pp_top;

  /// If true, adds d(Pc^up)/d(sat_n) * grad(sat_n) to the CO2 Darcy potential.
  /// Requires VEUpscaledCapPressure to be present (declares ve_dPcup_dsatn).
  const bool _capillary;

  /// d(Pc^up)/d(sat_n) [Pa] from VEUpscaledCapPressure; null when capillary=false.
  const MaterialProperty<Real> * const _dPcup_dsatn;
};
