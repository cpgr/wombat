#pragma once

#include "ADTimeKernelValue.h"
#include "VEDictator.h"

/**
 * VEMassTimeDerivative
 *
 * Depth-integrated mass storage term for one fluid phase:
 *
 *   ( H * phi_bar * rho_c * S_c_bar )_new  -  ( H * phi_bar * rho_c * S_c_bar )_old
 *   ------------------------------------------------------------------------------------
 *                                    dt
 *
 * where:
 *   H       = formation thickness        (ve_H, time-invariant)
 *   phi_bar = depth-averaged porosity    (ve_phi_bar, time-invariant)
 *   rho_c   = phase density              (ve_density, AD)
 *   S_c     = depth-averaged saturation  (ve_saturation, AD)
 *
 * Inherits from ADTimeKernelValue rather than ADKernelValue so that MOOSE
 * correctly identifies this as a time-derivative contribution. This matters
 * for predictor/corrector schemes and residual-norm reporting. The manual
 * (new - old) / dt differencing is used rather than _u_dot so that the full
 * nonlinear accumulation function H*phi*rho(p)*S(sat_n) is differenced
 * correctly -- important once VEFluidStateBrineCO2 introduces
 * pressure-dependent density.
 *
 * One instance per fluid phase in the input file. For the CO2 phase
 * (fluid_phase=0) assign to variable=sat_n; for the brine phase
 * (fluid_phase=1) assign to variable=pp_top.
 *
 * Jacobian assembled automatically via MOOSE AD -- no hand-coded derivatives.
 */
class VEMassTimeDerivative : public ADTimeKernelValue
{
public:
  static InputParameters validParams();
  VEMassTimeDerivative(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;

  /// VE configuration object.
  const VEDictator & _dictator;

  /// Fluid phase this kernel acts on (0 = CO2, 1 = brine).
  const unsigned int _fluid_phase;

  /// Formation thickness H [m] -- time-invariant, not AD.
  const MaterialProperty<Real> & _H;

  /// Depth-averaged porosity phi_bar [-] -- time-invariant, not AD.
  const MaterialProperty<Real> & _phi_bar;

  /// Depth-averaged phase densities [kg/m3] -- AD, size 2: [rho_CO2, rho_brine].
  const ADMaterialProperty<std::vector<Real>> & _density;

  /// Old depth-averaged phase densities [kg/m3] -- Real (previous time step).
  const MaterialProperty<std::vector<Real>> & _density_old;

  /// Depth-averaged phase saturations [-] -- AD, size 2: [S_n, S_w].
  const ADMaterialProperty<std::vector<Real>> & _saturation;

  /// Old depth-averaged phase saturations [-] -- Real (previous time step).
  const MaterialProperty<std::vector<Real>> & _saturation_old;
};
