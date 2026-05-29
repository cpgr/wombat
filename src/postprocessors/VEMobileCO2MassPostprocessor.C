#include "VEMobileCO2MassPostprocessor.h"

registerMooseObject("wombatApp", VEMobileCO2MassPostprocessor);

InputParameters
VEMobileCO2MassPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Integrates rho_co2 * phi_bar * H * (sat_n - sat_n_trap) over the domain "
      "to give the total mobile CO2 mass [kg].");

  params.addRequiredCoupledVar("sat_n",
                               "Depth-averaged CO2 saturation (primary variable).");
  params.addCoupledVar("sat_n_trap",
                       0.0,
                       "Depth-averaged residually trapped CO2 saturation [-].  "
                       "Couple to the AuxVariable written by VEHysteresis; "
                       "defaults to zero (all CO2 mobile) when not coupled.");
  params.addRequiredParam<Real>("rho_co2", "CO2 density [kg/m^3].");

  return params;
}

VEMobileCO2MassPostprocessor::VEMobileCO2MassPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _sat_n(coupledValue("sat_n")),
    _sat_n_trap(coupledValue("sat_n_trap")),
    _H(getMaterialProperty<Real>("ve_H")),
    _phi_bar(getMaterialProperty<Real>("ve_phi_bar")),
    _rho_co2(getParam<Real>("rho_co2"))
{
}

Real
VEMobileCO2MassPostprocessor::computeQpIntegral()
{
  const Real sat_mobile = std::max(0.0, _sat_n[_qp] - _sat_n_trap[_qp]);
  return _rho_co2 * _phi_bar[_qp] * _H[_qp] * sat_mobile;
}
