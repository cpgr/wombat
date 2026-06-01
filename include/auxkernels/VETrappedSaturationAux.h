#pragma once

#include "AuxKernel.h"

class VERelativePermeability;

/**
 * VETrappedSaturationAux
 *
 * Writes sat_n_trap, the depth-averaged residually trapped CO2 saturation, by
 * delegating to the relperm UserObject's trappedSaturation(sat_n, sat_n_max). The
 * result feeds the trapped / mobile CO2 mass postprocessors
 * (VETrappedCO2MassPostprocessor uses sat_n_trap; VEMobileCO2MassPostprocessor uses
 * sat_n - sat_n_trap), which already accept it.
 *
 * For a non-hysteretic UO trappedSaturation returns 0, so this aux is a no-op there.
 *
 * Couples to:
 *   - sat_n     : current depth-averaged CO2 saturation
 *   - sat_n_max : turning point from VESaturationMaxAux
 */
class VETrappedSaturationAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VETrappedSaturationAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  const VERelativePermeability & _relperm_uo;
  const VariableValue & _sat_n;
  const VariableValue & _sat_n_max;
};
