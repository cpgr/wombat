#include "ConstantFluidProperties.h"

registerMooseObject("wombatApp", ConstantFluidProperties);

InputParameters
ConstantFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addClassDescription(
      "Constant fluid density and viscosity, independent of pressure and temperature "
      "(zero derivatives). Verification-case stand-in for a full EOS, consumed by "
      "VEFluidProperties in place of CO2FluidProperties / BrineFluidProperties.");

  params.addRequiredRangeCheckedParam<Real>(
      "density", "density > 0", "Constant fluid density [kg/m3].");
  params.addRequiredRangeCheckedParam<Real>(
      "viscosity", "viscosity > 0", "Constant dynamic viscosity [Pa*s].");

  return params;
}

ConstantFluidProperties::ConstantFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _density(getParam<Real>("density")),
    _viscosity(getParam<Real>("viscosity"))
{
}

Real
ConstantFluidProperties::rho_from_p_T(Real /*p*/, Real /*T*/) const
{
  return _density;
}

void
ConstantFluidProperties::rho_from_p_T(
    Real /*p*/, Real /*T*/, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = _density;
  drho_dp = 0.0;
  drho_dT = 0.0;
}

Real
ConstantFluidProperties::mu_from_p_T(Real /*p*/, Real /*T*/) const
{
  return _viscosity;
}

void
ConstantFluidProperties::mu_from_p_T(
    Real /*p*/, Real /*T*/, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu = _viscosity;
  dmu_dp = 0.0;
  dmu_dT = 0.0;
}
