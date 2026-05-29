#include "VEPermeability.h"

registerMooseObject("wombatApp", VEPermeability);

InputParameters
VEPermeability::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Provides depth-integrated upscaled permeability tensor ve_K_up [m3]. "
      "Each component accepts a scalar constant or a spatially varying AuxVariable.");

  params.addRequiredCoupledVar("K_up_xx", "K_up tensor xx-component [m3].");
  params.addRequiredCoupledVar("K_up_yy", "K_up tensor yy-component [m3].");
  params.addCoupledVar("K_up_xy", 0.0, "K_up tensor off-diagonal xy-component [m3].");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPermeability::VEPermeability(const InputParameters & parameters)
  : Material(parameters),
    _K_up_xx(coupledValue("K_up_xx")),
    _K_up_xy(coupledValue("K_up_xy")),
    _K_up_yy(coupledValue("K_up_yy")),
    _K_up(declareProperty<RealTensorValue>("ve_K_up"))
{
}

void
VEPermeability::computeQpProperties()
{
  _K_up[_qp] = RealTensorValue(_K_up_xx[_qp], _K_up_xy[_qp], 0.0,
                                _K_up_xy[_qp], _K_up_yy[_qp], 0.0,
                                0.0,           0.0,            0.0);
}
