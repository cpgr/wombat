#include "VEFVFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("wombatApp", VEFVFluidProperties);

InputParameters
VEFVFluidProperties::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "FV adapter publishing phase density/viscosity as the functors ve_density_n/w and "
      "ve_viscosity_n/w for the FV flux kernel, from two SinglePhaseFluidProperties "
      "UserObjects (fp_nw / fp_w) evaluated at the pp_top functor and a constant "
      "temperature. Face-correct counterpart of VEFluidProperties for the FV advective "
      "flux; VEFluidProperties is still needed for the elemental FV mass-storage kernel.");

  params.addRequiredParam<UserObjectName>(
      "fp_nw",
      "SinglePhaseFluidProperties UserObject for the non-wetting phase (phase 0, "
      "typically CO2).");
  params.addRequiredParam<UserObjectName>(
      "fp_w",
      "SinglePhaseFluidProperties UserObject for the wetting phase (phase 1, "
      "typically brine).");
  params.addRequiredParam<MooseFunctorName>(
      "pp_top", "Pore pressure at the top of the formation [Pa]; the EOS evaluation pressure.");
  params.addRangeCheckedParam<Real>(
      "temperature",
      313.15,
      "temperature > 0",
      "Formation temperature [K], used only to evaluate the fluid properties (isothermal). "
      "Ignored by ConstantFluidProperties.");

  return params;
}

VEFVFluidProperties::VEFVFluidProperties(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _pp_top(getFunctor<ADReal>("pp_top")),
    _temperature(getParam<Real>("temperature")),
    _fp_nw(getUserObject<SinglePhaseFluidProperties>("fp_nw")),
    _fp_w(getUserObject<SinglePhaseFluidProperties>("fp_w"))
{
  // Temperature is constant, so it enters as an ADReal with zero derivatives; the AD
  // wrt pp_top comes from evaluating the pp_top functor at the supplied face argument.
  addFunctorProperty<ADReal>(
      "ve_density_n",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_nw.rho_from_p_T(_pp_top(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_density_w",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_w.rho_from_p_T(_pp_top(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_viscosity_n",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_nw.mu_from_p_T(_pp_top(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_viscosity_w",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_w.mu_from_p_T(_pp_top(r, t), ADReal(_temperature)); });
}
