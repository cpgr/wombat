#include "VEFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("wombatApp", VEFluidProperties);

InputParameters
VEFluidProperties::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Two-phase (non-wetting / wetting) fluid properties from FluidProperties "
      "UserObjects. Provides ve_density and ve_viscosity (size 2: [non-wetting, wetting]) "
      "evaluated at the top-surface pore pressure pp_top and a constant temperature. Use "
      "ConstantFluidProperties for verification cases and CO2FluidProperties / "
      "BrineFluidProperties (or any SinglePhaseFluidProperties) for a full EOS, without "
      "touching the kernels.");

  params.addRequiredParam<UserObjectName>(
      "fp_nw",
      "SinglePhaseFluidProperties UserObject for the non-wetting phase (phase 0, "
      "typically CO2).");
  params.addRequiredParam<UserObjectName>(
      "fp_w",
      "SinglePhaseFluidProperties UserObject for the wetting phase (phase 1, "
      "typically brine).");
  params.addRequiredCoupledVar(
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
      "top_surface (default): evaluate rho/mu at pp_top, i.e. at the formation top z_T "
      "-- reproduces existing golds bit-for-bit. interface: evaluate at the CO2-brine "
      "contact z_T + h, using the hydrostatic pressure p = pp_top + rho_n(pp_top)*|g|*h "
      "through the CO2 column, with the sharp-interface thickness h = sat_n*H/(1-S_wr) "
      "clamped to [0, H]. The increment uses the top-surface CO2 density (a single, "
      "non-iterated hydrostatic step; CO2 compressibility over a thin plume column is "
      "small). Requires sat_n and S_wr; reads ve_H from VEGeometry.");

  params.addCoupledVar(
      "sat_n",
      "Depth-averaged CO2 saturation. Required only for eos_reference_depth=interface, "
      "where it sets the sharp-interface plume thickness h = sat_n*H/(1-S_wr).");
  params.addCoupledVar(
      "z_top",
      "Top-surface elevation z_T [m]. Required only for eos_reference_depth=interface "
      "(with z_bottom, gives H = z_top - z_bottom). Coupled directly here -- not read "
      "from VEGeometry's ve_H -- so interface mode works in FV models too, which build "
      "H inline and do not run VEGeometry.");
  params.addCoupledVar(
      "z_bottom",
      "Bottom-surface elevation z_B [m]. Required only for eos_reference_depth=interface "
      "(with z_top, gives H = z_top - z_bottom).");
  params.addRangeCheckedParam<Real>(
      "S_wr",
      "S_wr >= 0 & S_wr < 1",
      "Residual water saturation in the CO2 zone [-]. Required only for "
      "eos_reference_depth=interface; must match VEPlumeReconstruction.");
  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only the magnitude |g| is used, for the hydrostatic "
      "pressure increment in eos_reference_depth=interface mode.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEFluidProperties::VEFluidProperties(const InputParameters & parameters)
  : Material(parameters),
    _eos_at_interface(getParam<MooseEnum>("eos_reference_depth") == "interface"),
    _pp_top(adCoupledValue("pp_top")),
    _temperature(getParam<Real>("temperature")),
    _sat_n((_eos_at_interface && isCoupled("sat_n")) ? &adCoupledValue("sat_n") : nullptr),
    _S_wr(isParamValid("S_wr") ? getParam<Real>("S_wr") : 0.0),
    _gravity_magnitude(getParam<RealVectorValue>("gravity").norm()),
    _z_top((_eos_at_interface && isCoupled("z_top")) ? &coupledValue("z_top") : nullptr),
    _z_bottom((_eos_at_interface && isCoupled("z_bottom")) ? &coupledValue("z_bottom") : nullptr),
    _fp_nw(getUserObject<SinglePhaseFluidProperties>("fp_nw")),
    _fp_w(getUserObject<SinglePhaseFluidProperties>("fp_w")),
    _density(declareADProperty<std::vector<Real>>("ve_density")),
    _viscosity(declareADProperty<std::vector<Real>>("ve_viscosity"))
{
  if (_eos_at_interface)
    for (const auto & p : {"sat_n", "z_top", "z_bottom"})
      if (!isCoupled(p))
        paramError(p, "Must be supplied when eos_reference_depth = interface.");
  if (_eos_at_interface && !isParamValid("S_wr"))
    paramError("S_wr", "Must be supplied when eos_reference_depth = interface.");
}

void
VEFluidProperties::fillQpValues(bool at_interface)
{
  // Constant temperature -> zero derivatives; promoted to ADReal so the EOS AD
  // overloads can be called uniformly.
  const ADReal T = _temperature;
  const ADReal p_top = _pp_top[_qp];

  ADReal p = p_top;
  if (at_interface)
  {
    // Sharp-interface plume thickness, clamped to [0, H]. Self-contained: H comes from
    // the coupled z_top/z_bottom (NOT VEGeometry's ve_H) so this works unchanged in FV
    // models, which build H inline and do not run VEGeometry; and h is built from sat_n
    // here rather than read from VEPlumeReconstruction's ve_h (which would form a
    // material-dependency cycle -- VEPlumeReconstruction already consumes ve_density).
    // The reference depth is itself a modeling approximation, so the sharp-interface h
    // is more than adequate: it differs from the capillary-fringe h only at second order
    // in the (small) CO2 compressibility.
    const Real H = (*_z_top)[_qp] - (*_z_bottom)[_qp];
    ADReal h = (*_sat_n)[_qp] * H / (1.0 - _S_wr);
    if (MetaPhysicL::raw_value(h) <= 0.0)
      h = 0.0;
    else if (MetaPhysicL::raw_value(h) >= H)
      h = H;

    // Hydrostatic pressure at the CO2-brine contact z_T + h, going down through the
    // CO2 column. rho_n is taken at the top-surface pressure (single non-iterated
    // step). p carries AD wrt pp_top (explicit term + rho_n_top) and wrt sat_n (h).
    const ADReal rho_n_top = _fp_nw.rho_from_p_T(p_top, T);
    p = p_top + rho_n_top * _gravity_magnitude * h;
  }

  _density[_qp].resize(2);
  _density[_qp][0] = _fp_nw.rho_from_p_T(p, T);
  _density[_qp][1] = _fp_w.rho_from_p_T(p, T);

  _viscosity[_qp].resize(2);
  _viscosity[_qp][0] = _fp_nw.mu_from_p_T(p, T);
  _viscosity[_qp][1] = _fp_w.mu_from_p_T(p, T);
}

void
VEFluidProperties::initQpStatefulProperties()
{
  // Seed old-value storage so the mass-storage kernels see correct densities at the
  // first timestep and do not compute a spurious initial storage spike. This uses the
  // same reference depth as computeQpProperties: H now comes from the coupled
  // z_top/z_bottom and h from sat_n -- all coupled VARIABLE values, which (unlike
  // material properties) are populated from the ICs before stateful init runs, so the
  // interface pressure is available here too.
  fillQpValues(_eos_at_interface);
}

void
VEFluidProperties::computeQpProperties()
{
  fillQpValues(_eos_at_interface);
}
