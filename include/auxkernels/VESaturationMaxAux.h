#pragma once

#include "AuxKernel.h"

/**
 * VESaturationMaxAux
 *
 * Stateful auxiliary that tracks sat_n_max, the maximum CO2 saturation each point
 * has ever experienced -- the turning point of the imbibition scanning curve used by
 * the hysteretic relative permeability (VERelPermHysteresisUO) for residual trapping.
 *
 *   sat_n_max <- max(sat_n, sat_n_max_old)
 *
 * advanced once per step (execute_on = TIMESTEP_END, after the nonlinear solve
 * converges). It is therefore LAGGED / explicit: within a solve sat_n_max is frozen,
 * which keeps the residual smooth and matches the "hysteresis is state, not a derived
 * quantity" design. The aux variable must be given an initial condition equal to the
 * sat_n IC, and declared with the same family as sat_n (LAGRANGE for FE, MONOMIAL /
 * MooseVariableFVReal for FV).
 *
 * Couples to:
 *   - sat_n : the depth-averaged CO2 saturation (primary variable)
 * and reads its own old value via uOld().
 */
class VESaturationMaxAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VESaturationMaxAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  /// Current depth-averaged CO2 saturation.
  const VariableValue & _sat_n;

  /// Previous-step value of sat_n_max (this aux's own old solution).
  const VariableValue & _sat_n_max_old;
};
