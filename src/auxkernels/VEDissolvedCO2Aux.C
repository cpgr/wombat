#include "VEDissolvedCO2Aux.h"

registerMooseObject("wombatApp", VEDissolvedCO2Aux);

InputParameters
VEDissolvedCO2Aux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Accumulates the areal dissolved CO2 mass c_diss (kg/m^2): "
      "c_diss <- c_diss_old + ve_dissolution_rate*dt, advanced at TIMESTEP_END with the "
      "converged rate so it matches the implicit VEDissolutionSink (mass conservation). "
      "Integrate over the domain for the total dissolved CO2 mass.");

  // Advance after the solve converges; lagged/explicit (frozen within the next solve).
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;

  return params;
}

VEDissolvedCO2Aux::VEDissolvedCO2Aux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _rate(getADMaterialProperty<Real>("ve_dissolution_rate")),
    _c_diss_old(uOld())
{
}

Real
VEDissolvedCO2Aux::computeValue()
{
  return _c_diss_old[_qp] + MetaPhysicL::raw_value(_rate[_qp]) * _dt;
}
