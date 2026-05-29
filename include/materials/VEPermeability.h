#pragma once

#include "Material.h"

/**
 * VEPermeability
 *
 * Provides the depth-integrated, upscaled permeability tensor ve_K_up [m3]
 * (units: m2 * m because H is already folded in by the upscaling workflow).
 *
 * Each independent component of the symmetric 2x2 tensor is supplied via a
 * required or optional coupled variable, accepting either a scalar constant or
 * a spatially varying AuxVariable from the Exodus mesh:
 *
 *   K_up_xx  (required)  -- xx-component [m3]
 *   K_up_yy  (required)  -- yy-component [m3]
 *   K_up_xy  (optional, default 0)  -- off-diagonal component [m3]
 *
 * The tensor is symmetric; K_up_zz = 0 (no vertical flow in a VE model).
 *
 * ve_K_up carries no Jacobian entries (permeability is time-invariant).
 */
class VEPermeability : public Material
{
public:
  static InputParameters validParams();
  VEPermeability(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const VariableValue & _K_up_xx;
  const VariableValue & _K_up_xy;
  const VariableValue & _K_up_yy;

  MaterialProperty<RealTensorValue> & _K_up;
};
