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

  /// If true, adds grad(Pc^up) = d(Pc^up)/d(sat_n)*grad(sat_n) +
  /// d(Pc^up)/d(H)*grad(H) to the CO2 Darcy potential. Requires
  /// VEUpscaledCapPressure (ve_dPcup_dsatn, ve_dPcup_dH) and VEGeometry (ve_grad_H).
  const bool _capillary;

  /// d(Pc^up)/d(sat_n) [Pa] from VEUpscaledCapPressure; null when capillary=false.
  /// AD: in sharp_interface mode it is constant in sat_n (zero derivative -> the
  /// grad(sat_n) Jacobian rides on _grad_u alone, as before); in capillary_fringe
  /// mode it depends on sat_n through S_n(h), and its AD derivative makes the
  /// grad(sat_n) term's Jacobian exact.
  const ADMaterialProperty<Real> * const _dPcup_dsatn;

  /// d(Pc^up)/d(H) [Pa/m] from VEUpscaledCapPressure; null when capillary=false.
  /// AD so the grad(H) term's sat_n Jacobian is captured (grad(H) is fixed geometry).
  const ADMaterialProperty<Real> * const _dPcup_dH;

  /// grad(H) [m/m] from VEGeometry (ve_grad_H); null when capillary=false.
  const MaterialProperty<RealGradient> * const _grad_H;
};
