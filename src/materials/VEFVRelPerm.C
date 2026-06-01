#include "VEFVRelPerm.h"
#include "VERelativePermeability.h"

registerMooseObject("wombatApp", VEFVRelPerm);

InputParameters
VEFVRelPerm::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Upscaled relative permeability as the functors ve_relperm_n / ve_relperm_w for the "
      "FV flux kernel, delegating the kr_c(sat_n) curve to a VERelativePermeability "
      "UserObject. Evaluated at the upwind face so Dirichlet saturation inlets drive "
      "inflow. FE counterpart: VERelPerm.");

  params.addRequiredParam<UserObjectName>(
      "relperm_uo", "VERelativePermeability UserObject providing the kr_c(sat_n) curve.");
  params.addRequiredParam<MooseFunctorName>(
      "sat_n", "Depth-averaged CO2 (non-wetting) saturation -- the FV primary variable.");
  params.addParam<MooseFunctorName>(
      "sat_n_max",
      "Optional turning-point functor (max historical CO2 saturation, from "
      "VESaturationMaxAux) for hysteretic relperm models. Omit for non-hysteretic "
      "(drainage-only) curves.");

  return params;
}

VEFVRelPerm::VEFVRelPerm(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _relperm_uo(getUserObject<VERelativePermeability>("relperm_uo")),
    _sat_n(getFunctor<ADReal>("sat_n")),
    _has_sat_n_max(isParamValid("sat_n_max")),
    _sat_n_max(_has_sat_n_max ? &getFunctor<ADReal>("sat_n_max") : nullptr)
{
  addFunctorProperty<ADReal>(
      "ve_relperm_n",
      [this](const auto & r, const auto & t) -> ADReal
      {
        if (_has_sat_n_max)
          // sat_n_max is a lagged aux -> frozen; pass its raw value so the AD seed
          // stays wrt sat_n only.
          return _relperm_uo.relativePermeabilityAD(
              _sat_n(r, t), MetaPhysicL::raw_value((*_sat_n_max)(r, t)), 0);
        return _relperm_uo.relativePermeabilityAD(_sat_n(r, t), 0);
      });

  addFunctorProperty<ADReal>(
      "ve_relperm_w",
      [this](const auto & r, const auto & t) -> ADReal
      {
        if (_has_sat_n_max)
          return _relperm_uo.relativePermeabilityAD(
              _sat_n(r, t), MetaPhysicL::raw_value((*_sat_n_max)(r, t)), 1);
        return _relperm_uo.relativePermeabilityAD(_sat_n(r, t), 1);
      });
}
