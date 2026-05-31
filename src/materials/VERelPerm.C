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

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VERelPerm::VERelPerm(const InputParameters & parameters)
  : Material(parameters),
    _relperm_uo(getUserObject<VERelativePermeability>("relperm_uo")),
    _saturation(getADMaterialProperty<std::vector<Real>>("ve_saturation")),
    _relperm(declareADProperty<std::vector<Real>>("ve_relperm"))
{
}

void
VERelPerm::computeQpProperties()
{
  _relperm[_qp].resize(2);

  const ADReal & sat_n = _saturation[_qp][0];
  _relperm[_qp][0] = _relperm_uo.relativePermeabilityAD(sat_n, 0);
  _relperm[_qp][1] = _relperm_uo.relativePermeabilityAD(sat_n, 1);
}
