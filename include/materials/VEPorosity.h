#pragma once

#include "Material.h"

/**
 * VEPorosity
 *
 * Provides the depth-averaged porosity ve_phi_bar [-].
 *
 * Two modes selected by the 'use_coupled_var' parameter:
 *
 *   Constant mode (default, for verification cases):
 *     ve_phi_bar = phi_bar   (uniform scalar parameter)
 *
 *   Coupled mode (for field cases with Exodus-supplied porosity):
 *     ve_phi_bar = phi_bar_var   (nodal AuxVariable from the mesh)
 *
 * ve_phi_bar is a plain Real because porosity is treated as time-invariant
 * and carries no Jacobian entries.
 */
class VEPorosity : public Material
{
public:
  static InputParameters validParams();
  VEPorosity(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const bool _use_coupled;

  /// Uniform porosity value used in constant mode.
  const Real _phi_bar_const;

  /// Nodal AuxVariable value used in coupled mode (null in constant mode).
  const VariableValue * _phi_bar_var;

  MaterialProperty<Real> & _phi_bar;
};
