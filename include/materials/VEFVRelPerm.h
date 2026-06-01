#pragma once

#include "FunctorMaterial.h"

class VERelativePermeability;

/**
 * VEFVRelPerm
 *
 * FV adapter: publishes the upscaled relative permeabilities as the functors
 * ve_relperm_n / ve_relperm_w for the FV flux kernel, delegating the kr_c(sat_n)
 * curve to a VERelativePermeability UserObject.
 *
 * The functors evaluate sat_n at whatever face argument the flux kernel supplies
 * (the upwind face), so a Dirichlet sat_n inlet drives inflow correctly. The FE
 * counterpart is VERelPerm.
 *
 * For hysteretic models, supply the optional sat_n_max functor (the turning point
 * from VESaturationMaxAux); it is evaluated at the same face and passed FROZEN (raw
 * value) to the UO's scanning-curve overload. When sat_n_max is not supplied the
 * non-hysteretic two-argument curve is used, so existing inputs are unchanged.
 */
class VEFVRelPerm : public FunctorMaterial
{
public:
  static InputParameters validParams();
  VEFVRelPerm(const InputParameters & parameters);

protected:
  const VERelativePermeability & _relperm_uo;
  const Moose::Functor<ADReal> & _sat_n;

  /// Whether the turning point sat_n_max is supplied (hysteresis active).
  const bool _has_sat_n_max;
  /// Turning point functor; only dereferenced when _has_sat_n_max.
  const Moose::Functor<ADReal> * const _sat_n_max;
};
