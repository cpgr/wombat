#pragma once

#include "FVElementalKernel.h"

/**
 * VEFVMassTimeDerivative
 *
 * Finite-volume depth-integrated mass storage term for one fluid phase:
 *
 *   ( H * phi_bar * rho_c * S_c )_new  -  ( H * phi_bar * rho_c * S_c )_old
 *   -------------------------------------------------------------------------
 *                                   dt
 *
 * FV analogue of VEMassTimeDerivative.  Formation thickness H is built inline
 * from the FV-variable cell values z_top - z_bottom (same reasoning as
 * VEFVAdvectiveFlux: FE-coupled geometry materials are not reinitialised on FV
 * faces, so geometry is sourced from FV variables for consistency between the
 * storage and flux terms).  Use fluid_phase = 0 on variable = sat_n and
 * fluid_phase = 1 on variable = pp_top.
 */
class VEFVMassTimeDerivative : public FVElementalKernel
{
public:
  static InputParameters validParams();
  VEFVMassTimeDerivative(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const unsigned int _fluid_phase;

  // Geometry from FV-variable cell values (static, non-AD).
  const VariableValue & _z_top;
  const VariableValue & _z_bottom;

  const MaterialProperty<Real> & _phi_bar;

  const ADMaterialProperty<std::vector<Real>> & _density;
  const MaterialProperty<std::vector<Real>> & _density_old;

  const ADMaterialProperty<std::vector<Real>> & _saturation;
  const MaterialProperty<std::vector<Real>> & _saturation_old;
};
