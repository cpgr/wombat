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

  return params;
}

VEFVRelPerm::VEFVRelPerm(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _relperm_uo(getUserObject<VERelativePermeability>("relperm_uo")),
    _sat_n(getFunctor<ADReal>("sat_n"))
{
  addFunctorProperty<ADReal>(
      "ve_relperm_n",
      [this](const auto & r, const auto & t) -> ADReal
      { return _relperm_uo.relativePermeabilityAD(_sat_n(r, t), 0); });

  addFunctorProperty<ADReal>(
      "ve_relperm_w",
      [this](const auto & r, const auto & t) -> ADReal
      { return _relperm_uo.relativePermeabilityAD(_sat_n(r, t), 1); });
}
