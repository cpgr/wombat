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
 *
 * For hysteretic models, couple the optional sat_n_max AuxVariable (the turning
 * point from VESaturationMaxAux); it is passed FROZEN to the UO's scanning-curve
 * overload (the AD seed stays wrt sat_n only). When sat_n_max is not coupled the
 * non-hysteretic two-argument curve is used, so existing inputs are unchanged.
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

  /// Whether the turning point sat_n_max is coupled (hysteresis active).
  const bool _has_sat_n_max;
  /// Turning point sat_n_max [-] (frozen, lagged aux); only valid if _has_sat_n_max.
  const VariableValue & _sat_n_max;

  ADMaterialProperty<std::vector<Real>> & _relperm;
};
