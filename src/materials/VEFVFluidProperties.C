#include "VEFVFluidProperties.h"

registerMooseObject("wombatApp", VEFVFluidProperties);

InputParameters
VEFVFluidProperties::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "FV adapter publishing phase density/viscosity as the functors ve_density_n/w and "
      "ve_viscosity_n/w for the FV flux kernel, from two SinglePhaseFluidProperties "
      "UserObjects (fp_nw / fp_w) evaluated at a reference pressure (top-surface pp_top or "
      "the z_T + h interface pressure, see eos_reference_depth) and a constant temperature. "
      "Face-correct counterpart of VEFluidProperties for the FV advective flux; "
      "VEFluidProperties is still needed for the elemental FV mass-storage kernel.");

  params.addRequiredParam<UserObjectName>(
      "fp_nw",
      "SinglePhaseFluidProperties UserObject for the non-wetting phase (phase 0, "
      "typically CO2).");
  params.addRequiredParam<UserObjectName>(
      "fp_w",
      "SinglePhaseFluidProperties UserObject for the wetting phase (phase 1, "
      "typically brine).");
  params.addRequiredParam<MooseFunctorName>(
      "pp_top", "Pore pressure at the top of the formation [Pa]; the top-surface EOS pressure.");
  params.addRangeCheckedParam<Real>(
      "temperature",
      313.15,
      "temperature > 0",
      "Formation temperature [K], used only to evaluate the fluid properties (isothermal). "
      "Ignored by ConstantFluidProperties.");

  MooseEnum eos_reference_depth("top_surface interface", "top_surface");
  params.addParam<MooseEnum>(
      "eos_reference_depth",
      eos_reference_depth,
      "Depth at which the EOS pressure is evaluated (CLAUDE.md key-subtlety 2). "
      "top_surface (default): evaluate rho/mu at pp_top (formation top z_T) -- reproduces "
      "existing golds. interface: evaluate at the CO2-brine contact z_T + h, using the "
      "hydrostatic pressure p = pp_top + rho_n(pp_top)*|g|*h with the sharp-interface "
      "thickness h = sat_n*(z_top - z_bottom)/(1 - S_wr) clamped to [0, H]. Requires the "
      "sat_n, z_top and z_bottom functors. Mirrors VEFluidProperties.");

  params.addParam<MooseFunctorName>(
      "sat_n",
      "Depth-averaged CO2 saturation functor. Required only for "
      "eos_reference_depth=interface.");
  params.addParam<MooseFunctorName>(
      "z_top", "Top-surface elevation z_T functor [m]. Required only for "
               "eos_reference_depth=interface (with z_bottom, gives H).");
  params.addParam<MooseFunctorName>(
      "z_bottom", "Bottom-surface elevation z_B functor [m]. Required only for "
                  "eos_reference_depth=interface (with z_top, gives H).");
  params.addRangeCheckedParam<Real>(
      "S_wr",
      0.0,
      "S_wr >= 0 & S_wr < 1",
      "Residual water saturation in the CO2 zone [-]. Used only in "
      "eos_reference_depth=interface; must match VEPlumeReconstruction.");
  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only |g| is used, for the hydrostatic increment in "
      "eos_reference_depth=interface mode.");

  return params;
}

VEFVFluidProperties::VEFVFluidProperties(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _eos_at_interface(getParam<MooseEnum>("eos_reference_depth") == "interface"),
    _pp_top(getFunctor<ADReal>("pp_top")),
    _temperature(getParam<Real>("temperature")),
    _S_wr(getParam<Real>("S_wr")),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _sat_n((_eos_at_interface && isParamValid("sat_n")) ? &getFunctor<ADReal>("sat_n") : nullptr),
    _z_top((_eos_at_interface && isParamValid("z_top")) ? &getFunctor<ADReal>("z_top") : nullptr),
    _z_bottom((_eos_at_interface && isParamValid("z_bottom")) ? &getFunctor<ADReal>("z_bottom")
                                                              : nullptr),
    _fp_nw(getUserObject<SinglePhaseFluidProperties>("fp_nw")),
    _fp_w(getUserObject<SinglePhaseFluidProperties>("fp_w"))
{
  if (_eos_at_interface)
    for (const auto & p : {"sat_n", "z_top", "z_bottom"})
      if (!isParamValid(p))
        paramError(p, "Must be supplied when eos_reference_depth = interface.");

  // Temperature is constant, so it enters as an ADReal with zero derivatives; the AD
  // wrt pp_top (and, in interface mode, wrt sat_n) comes from evaluating the functors
  // inside refPressure at the supplied face argument.
  addFunctorProperty<ADReal>(
      "ve_density_n",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_nw.rho_from_p_T(refPressure(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_density_w",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_w.rho_from_p_T(refPressure(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_viscosity_n",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_nw.mu_from_p_T(refPressure(r, t), ADReal(_temperature)); });
  addFunctorProperty<ADReal>(
      "ve_viscosity_w",
      [this](const auto & r, const auto & t) -> ADReal
      { return _fp_w.mu_from_p_T(refPressure(r, t), ADReal(_temperature)); });
}
