#include "VEFluidPropertiesConst.h"

registerMooseObject("wombatApp", VEFluidPropertiesConst);

InputParameters
VEFluidPropertiesConst::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Constant fluid properties for VE verification cases. Provides ve_density and "
      "ve_viscosity with user-supplied scalar values for CO2 and brine.");

  params.addRequiredParam<Real>("rho_co2",   "CO2 density [kg/m3].");
  params.addRequiredParam<Real>("rho_brine", "Brine density [kg/m3].");
  params.addRequiredParam<Real>("mu_co2",    "CO2 dynamic viscosity [Pa*s].");
  params.addRequiredParam<Real>("mu_brine",  "Brine dynamic viscosity [Pa*s].");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEFluidPropertiesConst::VEFluidPropertiesConst(const InputParameters & parameters)
  : Material(parameters),
    _rho_co2(getParam<Real>("rho_co2")),
    _rho_brine(getParam<Real>("rho_brine")),
    _mu_co2(getParam<Real>("mu_co2")),
    _mu_brine(getParam<Real>("mu_brine")),
    _density(declareADProperty<std::vector<Real>>("ve_density")),
    _viscosity(declareADProperty<std::vector<Real>>("ve_viscosity"))
{
}

void
VEFluidPropertiesConst::fillQpValues()
{
  _density[_qp].resize(2);
  _density[_qp][0] = _rho_co2;
  _density[_qp][1] = _rho_brine;

  _viscosity[_qp].resize(2);
  _viscosity[_qp][0] = _mu_co2;
  _viscosity[_qp][1] = _mu_brine;
}

void
VEFluidPropertiesConst::initQpStatefulProperties()
{
  // Seed old-value storage so VEMassTimeDerivative sees correct densities
  // at the first timestep and does not compute a spurious initial storage spike.
  fillQpValues();
}

void
VEFluidPropertiesConst::computeQpProperties()
{
  fillQpValues();
}
