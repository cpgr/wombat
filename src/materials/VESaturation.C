#include "VESaturation.h"

registerMooseObject("wombatApp", VESaturation);

InputParameters
VESaturation::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Wraps the primary variable sat_n into the two-element AD saturation vector "
      "ve_saturation = (sat_n, 1 - sat_n) consumed by VE kernels and materials.");

  params.addRequiredCoupledVar(
      "sat_n",
      "Depth-averaged CO2 (non-wetting) saturation -- primary variable of the simulation.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VESaturation::VESaturation(const InputParameters & parameters)
  : Material(parameters),
    _sat_n(adCoupledValue("sat_n")),
    _sat_n_real(coupledValue("sat_n")),
    _saturation(declareADProperty<std::vector<Real>>("ve_saturation"))
{
}

void
VESaturation::initQpStatefulProperties()
{
  // Use the plain Real value at the initial condition to avoid seeding
  // spurious AD derivatives into the old-value storage.
  _saturation[_qp].resize(2);
  _saturation[_qp][0] = _sat_n_real[_qp];
  _saturation[_qp][1] = 1.0 - _sat_n_real[_qp];
}

void
VESaturation::computeQpProperties()
{
  _saturation[_qp].resize(2);
  _saturation[_qp][0] = _sat_n[_qp];
  _saturation[_qp][1] = 1.0 - _sat_n[_qp];
}
