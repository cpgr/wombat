#include "VEPorosity.h"

registerMooseObject("wombatApp", VEPorosity);

InputParameters
VEPorosity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Provides depth-averaged porosity ve_phi_bar. "
      "Supply 'phi_bar' as a constant (e.g. phi_bar = 0.2) or couple it to "
      "a spatially varying AuxVariable from the Exodus mesh (field cases).");

  params.addRequiredCoupledVar("phi_bar", "Depth-averaged porosity (-).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPorosity::VEPorosity(const InputParameters & parameters)
  : Material(parameters),
    _phi_bar_var(coupledValue("phi_bar")),
    _phi_bar(declareProperty<Real>("ve_phi_bar"))
{
}

void
VEPorosity::computeQpProperties()
{
  _phi_bar[_qp] = _phi_bar_var[_qp];
}
