#pragma once

#include "AuxKernel.h"

/**
 * VEDissolvedCO2Aux
 *
 * Accumulates the areal dissolved CO2 mass c_diss [kg/m^2] consumed from the free phase by
 * VEDissolutionSink:
 *
 *   c_diss <- c_diss_old + ve_dissolution_rate * dt
 *
 * advanced once per step (execute_on = TIMESTEP_END) using the converged end-of-step rate,
 * so it matches the implicit (backward-Euler) sink exactly -> the dissolved inventory
 * gained equals the free CO2 mass removed (mass conservation). Lagged/explicit: the value
 * is frozen within the next solve, so it can feed the VEDissolution capacity cap without
 * adding a nonlinear coupling.
 *
 * The total dissolved CO2 mass is the domain integral of c_diss
 * (ElementIntegralVariablePostprocessor). Sink-only / Option A: c_diss is an immobile
 * inventory, not transported.
 *
 * The c_diss variable must be ELEMENTAL (CONSTANT MONOMIAL): this aux reads the
 * ve_dissolution_rate material property, which is defined at quadrature points, not nodes.
 */
class VEDissolvedCO2Aux : public AuxKernel
{
public:
  static InputParameters validParams();
  VEDissolvedCO2Aux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  /// Areal dissolution rate [kg/m^2/s] from VEDissolution.
  const ADMaterialProperty<Real> & _rate;
  /// Previous-step value of c_diss (this aux's own old solution).
  const VariableValue & _c_diss_old;
};
