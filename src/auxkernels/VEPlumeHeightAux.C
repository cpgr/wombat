#include "VEPlumeHeightAux.h"

registerMooseObject("wombatApp", VEPlumeHeightAux);

InputParameters
VEPlumeHeightAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes CO2 plume thickness h = sat_n * H / (1 - S_wr) under the "
      "sharp-interface VE assumption.  Result is clamped to [0, H].");

  params.addRequiredCoupledVar("sat_n",
                               "Depth-averaged CO2 saturation (primary variable).");

  params.addParam<Real>("S_wr",
                        0.0,
                        "Residual water saturation in the CO2 zone [-].  "
                        "Set to zero for the basic Buckley-Leverett case.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEPlumeHeightAux::VEPlumeHeightAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _sat_n(coupledValue("sat_n")),
    _H(getMaterialProperty<Real>("ve_H")),
    _S_wr(getParam<Real>("S_wr"))
{
  if (_S_wr < 0.0 || _S_wr >= 1.0)
    paramError("S_wr", "S_wr must be in [0, 1).");
}

Real
VEPlumeHeightAux::computeValue()
{
  const Real h = _sat_n[_qp] * _H[_qp] / (1.0 - _S_wr);
  return std::max(0.0, std::min(h, _H[_qp]));
}
