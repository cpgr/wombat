#pragma once

#include "Material.h"

/**
 * VEPermeability
 *
 * Provides the depth-integrated, upscaled permeability tensor ve_K_up [m3]
 * (units: m2 * m because H is already folded in by the upscaling workflow).
 *
 * Two modes:
 *
 *   Isotropic scalar mode (default, for verification cases):
 *     K_up_xx = K_up_yy = K_up   (scalar parameter, off-diagonal = 0)
 *
 *   Anisotropic mode (for field cases):
 *     K_up_xx, K_up_xy, K_up_yy supplied as parameters or AuxVariables.
 *     The tensor is symmetric 2x2 embedded in a 3x3 with K_up_zz = 0.
 *
 * ve_K_up is a plain Real tensor -- permeability carries no Jacobian entries
 * with respect to the primary variables (pp_top, sat_n).
 */
class VEPermeability : public Material
{
public:
  static InputParameters validParams();
  VEPermeability(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const bool _anisotropic;

  const Real _K_up_xx;
  const Real _K_up_xy;
  const Real _K_up_yy;

  MaterialProperty<RealTensorValue> & _K_up;
};
