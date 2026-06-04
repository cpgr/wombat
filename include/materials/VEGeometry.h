#pragma once

#include "Material.h"

/**
 * VEGeometry
 *
 * Computes two geometric quantities consumed by all other VE materials and kernels:
 *
 *   ve_H          -- formation thickness [m]
 *   ve_grad_z_top -- 2-D horizontal gradient of the top-surface elevation [m/m]
 *   ve_grad_H     -- 2-D horizontal gradient of the formation thickness [m/m]
 *
 * ve_grad_z_top enters every VEAdvectiveFlux kernel as the buoyancy drive
 * rho_c * g * grad(z_T).  ve_grad_H supplies the sat_n * grad(H) chain-rule
 * part of grad(Pc^up) for the capillary-diffusion term (consumed by
 * VEAdvectiveFluxS when capillary=true).  All three properties are plain Real
 * (not AD) because formation geometry is fixed and carries no Jacobian with
 * respect to the primary flow variables.
 *
 * Inputs: z_top and z_bottom (both nodal AuxVariables).  ve_H = z_top - z_bottom.
 * This is the single supported source of formation geometry; the ex2ve / Petrel
 * upscaling workflow emits nodal z_top / z_bottom directly, and coupling both
 * surfaces (rather than a precomputed thickness field) keeps geometry sourcing
 * identical in FE and FV and compatible with eos_reference_depth = interface,
 * which needs both surfaces.
 *
 * ve_grad_z_top is computed via coupledGradient("z_top"), which interpolates the
 * nodal z_top AuxVariable to quadrature points using the same LAGRANGE shape
 * functions.  z_top must therefore be declared as a LAGRANGE FIRST ORDER
 * AuxVariable; CONSTANT MONOMIAL would produce a zero gradient and silently kill
 * the entire buoyancy drive.  ve_grad_H is formed the same way as
 * grad(z_T) - grad(z_B), so z_bottom must likewise be LAGRANGE FIRST ORDER for a
 * non-zero grad(H); a zero grad(H) silently drops the topography part of the
 * capillary drive but leaves the constant-thickness physics exact.
 */
class VEGeometry : public Material
{
public:
  static InputParameters validParams();
  VEGeometry(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  // --- Inputs ---
  const VariableValue &    _z_top;
  const VariableGradient & _grad_z_top_var;
  const VariableValue &    _z_bottom;
  const VariableGradient & _grad_z_bottom_var;

  // --- Outputs ---
  MaterialProperty<Real> &         _H;
  MaterialProperty<RealGradient> & _grad_z_top;
  MaterialProperty<RealGradient> & _grad_H;
};
