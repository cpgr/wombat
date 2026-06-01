#include "VERelativePermeability.h"

InputParameters
VERelativePermeability::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Base class for VE upscaled relative-permeability models. Provides kr_c(sat_n) "
      "and its derivative as a swappable UserObject, consumed by the FE material "
      "VERelPerm and the FV functor material VEFVRelPerm.");
  return params;
}

VERelativePermeability::VERelativePermeability(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}

ADReal
VERelativePermeability::relativePermeabilityAD(const ADReal & sat_n, unsigned int phase) const
{
  const Real s = MetaPhysicL::raw_value(sat_n);
  ADReal kr = relativePermeability(s, phase);
  // Seed the sat_n derivative from the analytical d(kr)/d(sat_n) (chain rule).
  kr.derivatives() = dRelativePermeability(s, phase) * sat_n.derivatives();
  return kr;
}

ADReal
VERelativePermeability::relativePermeabilityAD(const ADReal & sat_n,
                                               Real sat_n_max,
                                               unsigned int phase) const
{
  const Real s = MetaPhysicL::raw_value(sat_n);
  ADReal kr = relativePermeability(s, sat_n_max, phase);
  // sat_n_max is frozen (lagged aux), so the only AD path is d(kr)/d(sat_n) along
  // the scanning curve through the fixed turning point sat_n_max.
  kr.derivatives() = dRelativePermeability(s, sat_n_max, phase) * sat_n.derivatives();
  return kr;
}
