#pragma once

#include "Material.h"

class SinglePhaseFluidProperties;

/**
 * VEFluidProperties
 *
 * Two-phase (non-wetting / wetting) fluid properties for the VE solver, evaluated
 * from FluidProperties-module UserObjects. The non-wetting phase is typically CO2
 * and the wetting phase brine, but nothing here assumes that -- any pair of
 * SinglePhaseFluidProperties works (hence the generic fp_nw / fp_w inputs, matching
 * the non-wetting/wetting convention of the relperm and capillary-pressure objects).
 *
 * Replaces the hard-coded VEFluidPropertiesConst: instead of four scalar inputs, the
 * density and viscosity of each phase come from a SinglePhaseFluidProperties
 * UserObject, so switching from constant properties to a full
 * pressure/temperature-dependent EOS is a UserObject swap in the input file -- the
 * same "plug and play" pattern PorousFlow uses (PorousFlowSingleComponentFluid +
 * SinglePhaseFluidProperties).
 *
 * Declares (identical shape to the old VEFluidPropertiesConst, so every kernel
 * and material that consumed those properties is unchanged), indexed by phase
 * (0 = non-wetting, 1 = wetting), matching the VEDictator phase convention:
 *   ve_density   [kg/m3]  -- ADMaterialProperty, size 2: [rho_nw, rho_w]
 *   ve_viscosity [Pa*s]   -- ADMaterialProperty, size 2: [mu_nw,  mu_w]
 *
 * Both are evaluated at the top-surface pore pressure pp_top and a constant
 * formation temperature. pp_top is coupled as an AD value, so for a genuine EOS
 * the d(rho)/d(pp_top) and d(mu)/d(pp_top) derivatives are seeded into the
 * Jacobian automatically; for ConstantFluidProperties they are zero, which is
 * correct and reproduces the verification golds bit-for-bit.
 *
 * Temperature is a constant scalar parameter (isothermal) used only to evaluate
 * the EOS; it is not a solution variable. The CLAUDE.md "EOS reference depth"
 * subtlety (evaluating density at the CO2-brine interface z_T + h rather than at
 * z_T) is a future switchable strategy; this material evaluates at pp_top.
 *
 * initQpStatefulProperties() seeds the old-value storage that
 * VEMassTimeDerivative / VEFVMassTimeDerivative read, so the storage term is
 * correct at t = 0.
 */
class VEFluidProperties : public Material
{
public:
  static InputParameters validParams();
  VEFluidProperties(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

private:
  void fillQpValues();

  /// Top-surface pore pressure [Pa] -- the EOS evaluation pressure (AD).
  const ADVariableValue & _pp_top;
  /// Constant formation temperature [K], used only to evaluate the EOS.
  const Real _temperature;

  const SinglePhaseFluidProperties & _fp_nw;
  const SinglePhaseFluidProperties & _fp_w;

  ADMaterialProperty<std::vector<Real>> & _density;
  ADMaterialProperty<std::vector<Real>> & _viscosity;
};
