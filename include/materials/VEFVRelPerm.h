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
 */
class VEFVRelPerm : public FunctorMaterial
{
public:
  static InputParameters validParams();
  VEFVRelPerm(const InputParameters & parameters);

protected:
  const VERelativePermeability & _relperm_uo;
  const Moose::Functor<ADReal> & _sat_n;
};
