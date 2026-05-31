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
