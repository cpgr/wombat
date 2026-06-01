#include "VERelPerm.h"
#include "VERelativePermeability.h"

registerMooseObject("wombatApp", VERelPerm);

InputParameters
VERelPerm::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Upscaled relative permeability as the qp material property ve_relperm for the FE "
      "kernels, delegating the kr_c(sat_n) curve to a VERelativePermeability UserObject. "
      "Reads ve_saturation. FV counterpart: VEFVRelPerm.");

  params.addRequiredParam<UserObjectName>(
      "relperm_uo", "VERelativePermeability UserObject providing the kr_c(sat_n) curve.");
  params.addCoupledVar(
      "sat_n_max",
      "Optional turning-point AuxVariable (max historical CO2 saturation, from "
      "VESaturationMaxAux) for hysteretic relperm models. Omit for non-hysteretic "
      "(drainage-only) curves.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VERelPerm::VERelPerm(const InputParameters & parameters)
  : Material(parameters),
    _relperm_uo(getUserObject<VERelativePermeability>("relperm_uo")),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _has_sat_n_max(isCoupled("sat_n_max")),
    _sat_n_max(_has_sat_n_max ? coupledValue("sat_n_max") : _zero),
    _relperm(declareADProperty<std::vector<Real>>("ve_relperm"))
{
}

void
VERelPerm::computeQpProperties()
{
  _relperm[_qp].resize(2);

  const ADReal & sat_n = _saturation[_qp][0];

  if (_has_sat_n_max)
  {
    // sat_n_max is a lagged aux -> frozen Real; AD seed stays wrt sat_n.
    const Real sat_n_max = _sat_n_max[_qp];
    _relperm[_qp][0] = _relperm_uo.relativePermeabilityAD(sat_n, sat_n_max, 0);
    _relperm[_qp][1] = _relperm_uo.relativePermeabilityAD(sat_n, sat_n_max, 1);
  }
  else
  {
    _relperm[_qp][0] = _relperm_uo.relativePermeabilityAD(sat_n, 0);
    _relperm[_qp][1] = _relperm_uo.relativePermeabilityAD(sat_n, 1);
  }
}
