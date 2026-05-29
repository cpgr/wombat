#pragma once

#include "Material.h"

/**
 * VEGeometry
 *
 * Reads the formation top and bottom surface elevations from coupled
 * AuxVariables (supplied by the Petrel-to-2D upscaling workflow via an
 * Exodus file) and computes the two geometric quantities that all other
 * VE materials and kernels depend on:
 *
 *   ve_H          = z_top - z_bottom   [m]
 *   ve_grad_z_top = grad(z_top)        [m/m]  (2D horizontal gradient)
 *
 * ve_grad_z_top drives the buoyancy migration term in VEAdvectiveFlux.
 * Both properties are plain Real (not AD) because formation geometry is
 * fixed -- it carries no Jacobian entries with respect to the primary
 * variables.
 */
class VEGeometry : public Material
{
public:
  static InputParameters validParams();
  VEGeometry(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  // --- Inputs: coupled AuxVariables from the Exodus mesh ---

  /// Top surface elevation z_top [m] (nodal AuxVar from Exodus).
  const VariableValue & _z_top;

  /// Gradient of z_top [m/m] -- MOOSE interpolates this from nodal values.
  const VariableGradient & _grad_z_top_var;

  /// Bottom surface elevation z_bottom [m] (nodal AuxVar from Exodus).
  const VariableValue & _z_bottom;

  // --- Outputs: material properties consumed by other VE objects ---

  MaterialProperty<Real> & _H;
  MaterialProperty<RealGradient> & _grad_z_top;
};
