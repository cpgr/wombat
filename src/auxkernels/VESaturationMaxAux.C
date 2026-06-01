#include "VESaturationMaxAux.h"

registerMooseObject("wombatApp", VESaturationMaxAux);

InputParameters
VESaturationMaxAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Stateful aux tracking sat_n_max = max(sat_n, sat_n_max_old), the turning point "
      "of the imbibition scanning curve for hysteretic residual trapping. Advanced once "
      "per step (lagged/explicit); requires an IC equal to the sat_n IC and the same "
      "variable family as sat_n.");

  params.addRequiredCoupledVar("sat_n",
                               "Depth-averaged CO2 saturation (the primary variable).");

  // Advance after the nonlinear solve converges so sat_n_max sees the converged sat_n
  // and stays frozen within the next solve.
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;

  return params;
}

VESaturationMaxAux::VESaturationMaxAux(const InputParameters & parameters)
  : AuxKernel(parameters), _sat_n(coupledValue("sat_n")), _sat_n_max_old(uOld())
{
}

Real
VESaturationMaxAux::computeValue()
{
  return std::max(_sat_n[_qp], _sat_n_max_old[_qp]);
}
