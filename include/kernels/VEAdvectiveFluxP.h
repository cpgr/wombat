#pragma once

#include "VEAdvectiveFlux.h"

/**
 * VEAdvectiveFluxP
 *
 * Depth-integrated Darcy advective flux kernel for use on the pressure
 * variable (variable = pp_top).
 *
 * Because the kernel variable IS pp_top, MOOSE's _grad_u[_qp] already
 * holds grad(pp_top) with full AD derivatives, so no extra coupling is
 * needed. This variant is used for the brine (wetting phase) mass
 * conservation equation:
 *
 *   fluid_phase = 1  (brine),  variable = pp_top
 */
class VEAdvectiveFluxP : public VEAdvectiveFluxBase
{
public:
  static InputParameters validParams();
  VEAdvectiveFluxP(const InputParameters & parameters);

protected:
  ADRealVectorValue precomputeQpResidual() override;
};
