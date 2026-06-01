#include "VETrappedSaturationAux.h"
#include "VERelativePermeability.h"

registerMooseObject("wombatApp", VETrappedSaturationAux);

InputParameters
VETrappedSaturationAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Residually trapped CO2 saturation sat_n_trap = "
      "relperm_uo.trappedSaturation(sat_n, sat_n_max), for the trapped / mobile CO2 "
      "mass postprocessors. Zero for a non-hysteretic relperm UO.");

  params.addRequiredParam<UserObjectName>(
      "relperm_uo",
      "VERelativePermeability UserObject providing trappedSaturation(sat_n, sat_n_max). "
      "Use VERelPermHysteresisUO for non-zero trapping.");
  params.addRequiredCoupledVar("sat_n", "Depth-averaged CO2 saturation (primary variable).");
  params.addRequiredCoupledVar(
      "sat_n_max", "Turning point (max historical CO2 saturation) from VESaturationMaxAux.");

  return params;
}

VETrappedSaturationAux::VETrappedSaturationAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _relperm_uo(getUserObject<VERelativePermeability>("relperm_uo")),
    _sat_n(coupledValue("sat_n")),
    _sat_n_max(coupledValue("sat_n_max"))
{
}

Real
VETrappedSaturationAux::computeValue()
{
  return _relperm_uo.trappedSaturation(_sat_n[_qp], _sat_n_max[_qp]);
}
