#pragma once

#include "Material.h"

class VERelativePermeability;

/**
 * VERelPerm
 *
 * FE adapter: publishes the upscaled relative permeabilities as the qp material
 * property ve_relperm (index 0 = CO2, 1 = brine) for the FE kernels, delegating
 * the kr_c(sat_n) curve to a VERelativePermeability UserObject.
 *
 * Reads ve_saturation (from VESaturation). The FV counterpart is VEFVRelPerm.
 */
class VERelPerm : public Material
{
public:
  static InputParameters validParams();
  VERelPerm(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const VERelativePermeability & _relperm_uo;
  const ADMaterialProperty<std::vector<Real>> & _saturation;
  ADMaterialProperty<std::vector<Real>> & _relperm;
};
