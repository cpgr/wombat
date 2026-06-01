#include "VEGeometry.h"

registerMooseObject("wombatApp", VEGeometry);

InputParameters
VEGeometry::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes ve_H (formation thickness) and ve_grad_z_top (top-surface gradient) "
      "from coupled AuxVariables. Two modes: top_bottom (z_top + z_bottom, computes H) "
      "and thickness (z_top + H supplied directly, e.g. from an Exodus upscaling mesh).");

  params.addRequiredCoupledVar(
      "z_top",
      "AuxVariable holding the top-surface elevation z_T [m]. "
      "Must be LAGRANGE FIRST ORDER so that coupledGradient produces a non-zero "
      "ve_grad_z_top for the buoyancy drive.");

  // top_bottom mode
  params.addCoupledVar(
      "z_bottom",
      "AuxVariable holding the bottom-surface elevation z_B [m]. "
      "Required in top_bottom mode; omit in thickness mode. "
      "Must be LAGRANGE FIRST ORDER for a non-zero grad(H) (see ve_grad_H).");

  // thickness mode
  params.addCoupledVar(
      "H",
      "AuxVariable holding the column-averaged formation thickness H [m]. "
      "Required in thickness mode; omit in top_bottom mode. "
      "Typical source: a nodal field on a Petrel-derived Exodus mesh. "
      "Must be LAGRANGE FIRST ORDER for a non-zero grad(H) (see ve_grad_H).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEGeometry::VEGeometry(const InputParameters & parameters)
  : Material(parameters),
    _thickness_mode(isCoupled("H")),
    _z_top(coupledValue("z_top")),
    _grad_z_top_var(coupledGradient("z_top")),
    _z_bottom(_thickness_mode ? nullptr : &coupledValue("z_bottom")),
    _grad_z_bottom_var(_thickness_mode ? nullptr : &coupledGradient("z_bottom")),
    _H_var(_thickness_mode ? &coupledValue("H") : nullptr),
    _grad_H_var(_thickness_mode ? &coupledGradient("H") : nullptr),
    _H(declareProperty<Real>("ve_H")),
    _grad_z_top(declareProperty<RealGradient>("ve_grad_z_top")),
    _grad_H(declareProperty<RealGradient>("ve_grad_H"))
{
  if (isCoupled("z_bottom") && isCoupled("H"))
    paramError("H",
               "Provide either z_bottom (top_bottom mode) or H (thickness mode), not both.");

  if (!isCoupled("z_bottom") && !isCoupled("H"))
    paramError("z_top",
               "Must provide either z_bottom (top_bottom mode) or H (thickness mode).");
}

void
VEGeometry::computeQpProperties()
{
  _H[_qp] = _thickness_mode ? (*_H_var)[_qp]
                             : _z_top[_qp] - (*_z_bottom)[_qp];

  _grad_z_top[_qp] = _grad_z_top_var[_qp];

  // grad(H) is the lateral thickness gradient. It feeds the sat_n*grad(H)
  // chain-rule part of grad(Pc^up) in the capillary-diffusion term (see
  // VEUpscaledCapPressure / VEFVCapPressure). thickness mode reads grad(H)
  // directly; top_bottom mode forms it as grad(z_T) - grad(z_B).
  if (_thickness_mode)
    _grad_H[_qp] = (*_grad_H_var)[_qp];
  else
    _grad_H[_qp] = _grad_z_top_var[_qp] - (*_grad_z_bottom_var)[_qp];
}
