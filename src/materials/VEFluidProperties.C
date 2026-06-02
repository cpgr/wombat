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
      "pp_top", "Pore pressure at the top of the formation [Pa]; the EOS evaluation pressure.");
  params.addRangeCheckedParam<Real>(
      "temperature",
      313.15,
      "temperature > 0",
      "Formation temperature [K], used only to evaluate the fluid properties (isothermal). "
      "Ignored by ConstantFluidProperties.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEFluidProperties::VEFluidProperties(const InputParameters & parameters)
  : Material(parameters),
    _pp_top(adCoupledValue("pp_top")),
    _temperature(getParam<Real>("temperature")),
    _fp_nw(getUserObject<SinglePhaseFluidProperties>("fp_nw")),
    _fp_w(getUserObject<SinglePhaseFluidProperties>("fp_w")),
    _density(declareADProperty<std::vector<Real>>("ve_density")),
    _viscosity(declareADProperty<std::vector<Real>>("ve_viscosity"))
{
}

void
VEFluidProperties::fillQpValues()
{
  const ADReal p = _pp_top[_qp];
  // Constant temperature -> zero derivatives; promoted to ADReal so the EOS AD
  // overloads can be called uniformly.
  const ADReal T = _temperature;

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
  // Seed old-value storage so the mass-storage kernels see correct densities at
  // the first timestep and do not compute a spurious initial storage spike.
  fillQpValues();
}

void
VEFluidProperties::computeQpProperties()
{
  fillQpValues();
}
