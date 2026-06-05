#include "VETrappedCO2MassPostprocessor.h"

registerMooseObject("wombatApp", VETrappedCO2MassPostprocessor);

InputParameters
VETrappedCO2MassPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Integrates rho_co2 * phi_bar * H * sat_n_trap over the domain "
      "to give the total residually trapped CO2 mass (kg).");

  params.addCoupledVar("sat_n_trap",
                       0.0,
                       "Depth-averaged residually trapped CO2 saturation (-).  "
                       "Couple to the AuxVariable written by VEHysteresis; "
                       "defaults to zero until hysteresis is active.");
  params.addRequiredParam<Real>("rho_co2", "CO2 density (kg/m^3).");

  return params;
}

VETrappedCO2MassPostprocessor::VETrappedCO2MassPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _sat_n_trap(coupledValue("sat_n_trap")),
    _H(getMaterialProperty<Real>("ve_H")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _rho_co2(getParam<Real>("rho_co2"))
{
}

Real
VETrappedCO2MassPostprocessor::computeQpIntegral()
{
  return _rho_co2 * _phi_bar[_qp] * _H[_qp] * std::max(0.0, _sat_n_trap[_qp]);
}
