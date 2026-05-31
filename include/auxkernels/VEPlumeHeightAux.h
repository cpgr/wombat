#pragma once

#include "AuxKernel.h"

/**
 * VEPlumeHeightAux
 *
 * Exposes the ve_h material property computed by VEPlumeReconstruction as an
 * elemental AuxVariable for output and postprocessing.
 *
 * Reads:
 *   - ve_h : material property from VEPlumeReconstruction (plume thickness [m])
 */
class VEPlumeHeightAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VEPlumeHeightAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  const ADMaterialProperty<Real> & _h;
};
