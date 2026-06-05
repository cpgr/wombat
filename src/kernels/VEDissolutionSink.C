#include "VEDissolutionSink.h"

registerMooseObject("wombatApp", VEDissolutionSink);

InputParameters
VEDissolutionSink::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription(
      "Sink on the CO2 mass equation removing free-phase CO2 at the convective dissolution "
      "rate ve_dissolution_rate (kg/m^2/s) (from VEDissolution). Assign to variable=sat_n. "
      "Residual contribution is +ve_dissolution_rate (a sink: free CO2 decreases). Off by "
      "default -- add only when dissolution is active.");

  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");

  return params;
}

VEDissolutionSink::VEDissolutionSink(const InputParameters & parameters)
  : ADKernelValue(parameters), _rate(getADMaterialProperty<Real>("ve_dissolution_rate"))
{
}

ADReal
VEDissolutionSink::precomputeQpResidual()
{
  // +rate: with storage d/dt(H phi rho_c S_c) + ... + rate = 0, the free CO2 mass falls.
  return _rate[_qp];
}
