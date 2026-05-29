#include "VEPorosity.h"

registerMooseObject("wombatApp", VEPorosity);

InputParameters
VEPorosity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Provides depth-averaged porosity ve_phi_bar. "
      "Use 'phi_bar' for a uniform scalar (verification) or couple 'phi_bar_var' "
      "for a spatially varying AuxVariable from the Exodus mesh (field cases).");

  params.addParam<Real>("phi_bar", 0.2, "Uniform depth-averaged porosity [-] (constant mode).");

  params.addParam<bool>("use_coupled_var",
                        false,
                        "If true, read porosity from the coupled AuxVariable 'phi_bar_var' "
                        "instead of the 'phi_bar' parameter.");

  params.addCoupledVar("phi_bar_var",
                       "AuxVariable holding depth-averaged porosity (coupled mode only).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPorosity::VEPorosity(const InputParameters & parameters)
  : Material(parameters),
    _use_coupled(getParam<bool>("use_coupled_var")),
    _phi_bar_const(getParam<Real>("phi_bar")),
    _phi_bar_var(_use_coupled ? &coupledValue("phi_bar_var") : nullptr),
    _phi_bar(declareProperty<Real>("ve_phi_bar"))
{
  if (_use_coupled && !isCoupled("phi_bar_var"))
    paramError("phi_bar_var",
               "use_coupled_var=true but no variable was coupled to phi_bar_var.");
}

void
VEPorosity::computeQpProperties()
{
  _phi_bar[_qp] = _use_coupled ? (*_phi_bar_var)[_qp] : _phi_bar_const;
}
