#pragma once

#include "Material.h"

/**
 * VEGeometry
 *
 * Computes two geometric quantities consumed by all other VE materials and kernels:
 *
 *   ve_H          -- formation thickness [m]
 *   ve_grad_z_top -- 2-D horizontal gradient of the top-surface elevation [m/m]
 *
 * ve_grad_z_top enters every VEAdvectiveFlux kernel as the buoyancy drive
 * rho_c * g * grad(z_T).  Both properties are plain Real (not AD) because
 * formation geometry is fixed and carries no Jacobian with respect to the
 * primary flow variables.
 *
 * Two input modes are supported; exactly one must be chosen:
 *
 *   top_bottom (default)
 *     Requires: z_top, z_bottom  (both nodal AuxVariables)
 *     Computes: ve_H = z_top - z_bottom
 *     Use for: synthetic / CSV-loaded geometries where the bottom surface
 *     elevation is explicitly available.
 *
 *   thickness
 *     Requires: z_top, H  (both nodal AuxVariables)
 *     Computes: ve_H = H  (used directly)
 *     Use for: real Exodus meshes from the Petrel upscaling workflow, where
 *     the workflow supplies H (column-averaged thickness) as a native field
 *     alongside phi_bar, K_up, etc.  z_bottom is not needed and need not be
 *     present on the mesh.
 *
 * In both modes ve_grad_z_top is computed via coupledGradient("z_top"),
 * which interpolates the nodal z_top AuxVariable to quadrature points using
 * the same LAGRANGE shape functions.  z_top must therefore be declared as a
 * LAGRANGE FIRST ORDER AuxVariable; CONSTANT MONOMIAL would produce a zero
 * gradient and silently kill the entire buoyancy drive.
 */
class VEGeometry : public Material
{
public:
  static InputParameters validParams();
  VEGeometry(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  // --- Input mode ---
  const bool _thickness_mode; ///< true = thickness mode, false = top_bottom mode

  // --- Inputs always required ---
  const VariableValue &    _z_top;
  const VariableGradient & _grad_z_top_var;

  // --- Inputs required in top_bottom mode ---
  const VariableValue * const _z_bottom; ///< nullptr in thickness mode

  // --- Inputs required in thickness mode ---
  const VariableValue * const _H_var; ///< nullptr in top_bottom mode

  // --- Outputs ---
  MaterialProperty<Real> &         _H;
  MaterialProperty<RealGradient> & _grad_z_top;
};
