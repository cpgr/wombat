#pragma once

#include "ADKernelGrad.h"

/**
 * VEAdvectiveFluxBase
 *
 * Abstract base class for the depth-integrated Darcy advective flux kernels.
 * After integration by parts the weak-form domain contribution for phase c is:
 *
 *   int grad(psi) . F_c dOmega
 *
 * where:
 *
 *   F_c = H * K_up * (kr_c_up / mu_c) * rho_c * (grad(pp_top) + rho_c * |g| * grad(z_T))
 *
 * The term rho_c * |g| * grad(z_T) is the sloped-top-surface buoyancy drive --
 * the dominant updip migration mechanism at field scale and a first-class
 * Jacobian term (carried through rho_c which is AD).
 *
 * This class is abstract because the pressure gradient grad(pp_top) is sourced
 * differently depending on which primary variable the kernel acts on:
 *
 *   VEAdvectiveFluxP  (variable = pp_top)
 *     grad(pp_top) = _grad_u[_qp]  (kernel's own gradient -- no extra coupling)
 *
 *   VEAdvectiveFluxS  (variable = sat_n)
 *     grad(pp_top) = adCoupledGradient("pp_top")[_qp]  (explicitly coupled)
 *
 * All material properties and parameters live here; derived classes only
 * implement precomputeQpResidual() with the correct pressure gradient source.
 */
class VEAdvectiveFluxBase : public ADKernelGrad
{
public:
  static InputParameters validParams();
  VEAdvectiveFluxBase(const InputParameters & parameters);

protected:
  /// Fluid phase this kernel acts on (0 = CO2, 1 = brine).
  const unsigned int _fluid_phase;

  /// Scalar magnitude of gravity [m/s2].
  const Real _gravity_magnitude;

  // --- Material properties ---

  const MaterialProperty<Real> & _H;
  const MaterialProperty<RealTensorValue> & _K_up;
  const ADMaterialProperty<std::vector<Real>> & _relperm;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _viscosity;
  const MaterialProperty<RealGradient> & _grad_z_top;
};
