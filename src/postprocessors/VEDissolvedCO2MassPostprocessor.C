#include "VEDissolvedCO2MassPostprocessor.h"

registerMooseObject("wombatApp", VEDissolvedCO2MassPostprocessor);

InputParameters
VEDissolvedCO2MassPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Integrates rho_brine * phi_bar * H * (1 - sat_n) * X_co2 over the domain "
      "to give the total dissolved CO2 mass (kg).");

  params.addRequiredCoupledVar("sat_n",
                               "Depth-averaged CO2 saturation (primary variable).");
  params.addCoupledVar("X_co2",
                       0.0,
                       "Dissolved CO2 mass fraction in brine (kg CO2 / kg brine).  "
                       "Couple to the AuxVariable written by VEDissolution; "
                       "defaults to zero until dissolution is active.");
  params.addRequiredParam<Real>("rho_brine", "Brine density (kg/m^3).");

  return params;
}

VEDissolvedCO2MassPostprocessor::VEDissolvedCO2MassPostprocessor(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _sat_n(coupledValue("sat_n")),
    _X_co2(coupledValue("X_co2")),
    _H(getMaterialProperty<Real>("ve_H")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _rho_brine(getParam<Real>("rho_brine"))
{
}

Real
VEDissolvedCO2MassPostprocessor::computeQpIntegral()
{
  const Real brine_sat = std::max(0.0, 1.0 - _sat_n[_qp]);
  return _rho_brine * _phi_bar[_qp] * _H[_qp] * brine_sat * std::max(0.0, _X_co2[_qp]);
}
