#include "VEGeometry.h"

registerMooseObject("wombatApp", VEGeometry);

InputParameters
VEGeometry::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes ve_H (formation thickness), ve_grad_z_top (top-surface gradient), and "
      "ve_grad_H (thickness gradient) from the coupled z_top / z_bottom AuxVariables: "
      "ve_H = z_top - z_bottom.");

  params.addRequiredCoupledVar(
      "z_top",
      "AuxVariable holding the top-surface elevation z_T (m). "
      "Must be LAGRANGE FIRST ORDER so that coupledGradient produces a non-zero "
      "ve_grad_z_top for the buoyancy drive.");

  params.addRequiredCoupledVar(
      "z_bottom",
      "AuxVariable holding the bottom-surface elevation z_B (m). "
      "Must be LAGRANGE FIRST ORDER for a non-zero grad(H) (see ve_grad_H).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEGeometry::VEGeometry(const InputParameters & parameters)
  : Material(parameters),
    _z_top(coupledValue("z_top")),
    _grad_z_top_var(coupledGradient("z_top")),
    _z_bottom(coupledValue("z_bottom")),
    _grad_z_bottom_var(coupledGradient("z_bottom")),
    _H(declareProperty<Real>("ve_H")),
    _grad_z_top(declareProperty<RealGradient>("ve_grad_z_top")),
    _grad_H(declareProperty<RealGradient>("ve_grad_H"))
{
}

void
VEGeometry::computeQpProperties()
{
  _H[_qp] = _z_top[_qp] - _z_bottom[_qp];

  _grad_z_top[_qp] = _grad_z_top_var[_qp];

  // grad(H) is the lateral thickness gradient, formed as grad(z_T) - grad(z_B).
  // It feeds the sat_n*grad(H) chain-rule part of grad(Pc^up) in the
  // capillary-diffusion term (see VEUpscaledCapPressure / VEFVCapPressure).
  _grad_H[_qp] = _grad_z_top_var[_qp] - _grad_z_bottom_var[_qp];
}
