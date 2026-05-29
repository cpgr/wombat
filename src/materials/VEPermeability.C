#include "VEPermeability.h"

registerMooseObject("wombatApp", VEPermeability);

InputParameters
VEPermeability::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Provides depth-integrated upscaled permeability tensor ve_K_up [m3]. "
      "Isotropic mode: set 'K_up' scalar. "
      "Anisotropic mode: set anisotropic=true and supply K_up_xx, K_up_xy, K_up_yy.");

  params.addParam<bool>("anisotropic",
                        false,
                        "If false (default), K_up_xx = K_up_yy = K_up and K_up_xy = 0. "
                        "If true, specify K_up_xx, K_up_xy, K_up_yy individually.");

  params.addParam<Real>("K_up",
                        1.0e-12,
                        "Isotropic depth-integrated permeability K_up [m3] (isotropic mode).");

  params.addParam<Real>("K_up_xx", 1.0e-12, "K_up tensor xx-component [m3] (anisotropic mode).");
  params.addParam<Real>("K_up_xy", 0.0, "K_up tensor xy-component [m3] (anisotropic mode).");
  params.addParam<Real>("K_up_yy", 1.0e-12, "K_up tensor yy-component [m3] (anisotropic mode).");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPermeability::VEPermeability(const InputParameters & parameters)
  : Material(parameters),
    _anisotropic(getParam<bool>("anisotropic")),
    _K_up_xx(_anisotropic ? getParam<Real>("K_up_xx") : getParam<Real>("K_up")),
    _K_up_xy(_anisotropic ? getParam<Real>("K_up_xy") : 0.0),
    _K_up_yy(_anisotropic ? getParam<Real>("K_up_yy") : getParam<Real>("K_up")),
    _K_up(declareProperty<RealTensorValue>("ve_K_up"))
{
}

void
VEPermeability::computeQpProperties()
{
  // Fill the symmetric 2x2 sub-block of a 3x3 tensor.
  // K_up_zz = 0 since there is no vertical flow in a VE model.
  _K_up[_qp] = RealTensorValue(_K_up_xx, _K_up_xy, 0.0,
                                _K_up_xy, _K_up_yy, 0.0,
                                0.0,      0.0,      0.0);
}
