#pragma once

#include "Material.h"

/**
 * VEPorosity
 *
 * Provides the depth-averaged porosity ve_phi_bar [-].
 *
 * Supply 'phi_bar' as a constant scalar (e.g. phi_bar = 0.2 for verification cases)
 * or couple it to a spatially varying AuxVariable from the Exodus mesh (field cases);
 * MOOSE resolves both through the same coupledValue interface. The parameter is required.
 *
 * ve_phi_bar carries no Jacobian entries (porosity is time-invariant).
 */
class VEPorosity : public Material
{
public:
  static InputParameters validParams();
  VEPorosity(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const VariableValue & _phi_bar_var;

  MaterialProperty<Real> & _phi_bar;
};
