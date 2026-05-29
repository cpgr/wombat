#include "VEGeometry.h"

registerMooseObject("wombatApp", VEGeometry);

InputParameters
VEGeometry::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes VE geometry properties ve_H and ve_grad_z_top from coupled "
      "AuxVariables z_top and z_bottom supplied by the Petrel upscaling workflow.");

  params.addRequiredCoupledVar("z_top",
                               "AuxVariable holding the top surface elevation z_T [m].");
  params.addRequiredCoupledVar("z_bottom",
                               "AuxVariable holding the bottom surface elevation z_B [m].");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEGeometry::VEGeometry(const InputParameters & parameters)
  : Material(parameters),
    _z_top(coupledValue("z_top")),
    _grad_z_top_var(coupledGradient("z_top")),
    _z_bottom(coupledValue("z_bottom")),
    _H(declareProperty<Real>("ve_H")),
    _grad_z_top(declareProperty<RealGradient>("ve_grad_z_top"))
{
}

void
VEGeometry::computeQpProperties()
{
  _H[_qp] = _z_top[_qp] - _z_bottom[_qp];
  _grad_z_top[_qp] = _grad_z_top_var[_qp];
}
