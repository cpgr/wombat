#include "VEFVDissolutionSink.h"

registerMooseObject("wombatApp", VEFVDissolutionSink);

InputParameters
VEFVDissolutionSink::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "FV sink on the CO2 mass equation removing free-phase CO2 at the convective "
      "dissolution rate ve_dissolution_rate (kg/m^2/s) (from VEDissolution). Assign to "
      "variable=sat_n. Residual contribution is +ve_dissolution_rate (a sink: free CO2 "
      "decreases). FV analogue of VEDissolutionSink; off by default -- add only when "
      "dissolution is active.");
  return params;
}

VEFVDissolutionSink::VEFVDissolutionSink(const InputParameters & parameters)
  : FVElementalKernel(parameters), _rate(getADMaterialProperty<Real>("ve_dissolution_rate"))
{
}

ADReal
VEFVDissolutionSink::computeQpResidual()
{
  // +rate: same elemental sign as VEFVMassTimeDerivative, so the free CO2 mass falls.
  return _rate[_qp];
}
