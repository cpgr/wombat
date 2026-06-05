#include "SimpleBrineFluidProperties.h"
#include "BrineFluidProperties.h"

registerMooseObject("wombatApp", SimpleBrineFluidProperties);

InputParameters
SimpleBrineFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addClassDescription(
      "Fixed-salinity brine as a two-variable (pressure, temperature) "
      "SinglePhaseFluidProperties. Wraps the FluidProperties-module BrineFluidProperties "
      "(pressure, temperature, salt mass fraction) at a constant salt_mass_fraction, "
      "exposing density and viscosity with derivatives wrt pressure and temperature only. "
      "Drop-in fp_w for VEFluidProperties when the wetting phase is brine rather than pure "
      "water.");

  params.addRequiredRangeCheckedParam<Real>(
      "salt_mass_fraction",
      "salt_mass_fraction >= 0 & salt_mass_fraction < 1",
      "Constant salt (NaCl) mass fraction of the brine (kg/kg).");
  params.addParam<UserObjectName>(
      "water_fp",
      "Optional water SinglePhaseFluidProperties UserObject passed through to the "
      "underlying BrineFluidProperties (e.g. a tabulated Water97). If omitted, "
      "BrineFluidProperties constructs its own water EOS.");

  return params;
}

SimpleBrineFluidProperties::SimpleBrineFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters), _xnacl(getParam<Real>("salt_mass_fraction"))
{
  // Build the underlying three-variable BrineFluidProperties UserObject, mirroring
  // how BrineFluidProperties itself constructs its component (water/NaCl) UOs.
  const std::string brine_name = name() + ":brine";
  {
    const std::string class_name = "BrineFluidProperties";
    InputParameters brine_params = _app.getFactory().getValidParams(class_name);
    if (parameters.isParamSetByUser("water_fp"))
      brine_params.set<UserObjectName>("water_fp") = getParam<UserObjectName>("water_fp");
    if (_tid == 0)
      _fe_problem.addUserObject(class_name, brine_name, brine_params);
  }
  _brine_fp = &_fe_problem.getUserObject<BrineFluidProperties>(brine_name);
}

Real
SimpleBrineFluidProperties::rho_from_p_T(Real p, Real T) const
{
  return _brine_fp->rho_from_p_T_X(p, T, _xnacl);
}

void
SimpleBrineFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  // Salinity is fixed, so the salinity derivative (drho_dx) is discarded.
  Real drho_dx;
  _brine_fp->rho_from_p_T_X(p, T, _xnacl, rho, drho_dp, drho_dT, drho_dx);
}

Real
SimpleBrineFluidProperties::mu_from_p_T(Real p, Real T) const
{
  return _brine_fp->mu_from_p_T_X(p, T, _xnacl);
}

void
SimpleBrineFluidProperties::mu_from_p_T(
    Real p, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  // Salinity is fixed, so the salinity derivative (dmu_dx) is discarded.
  Real dmu_dx;
  _brine_fp->mu_from_p_T_X(p, T, _xnacl, mu, dmu_dp, dmu_dT, dmu_dx);
}
