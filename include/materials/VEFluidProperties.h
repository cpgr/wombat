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
 * (0 = non-wetting, 1 = wetting), matching the VE phase convention:
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
 * the EOS; it is not a solution variable.
 *
 * EOS reference depth (CLAUDE.md key-subtlety 2), switchable via
 * eos_reference_depth:
 *   top_surface (default) -- evaluate rho/mu at pp_top (the formation top z_T).
 *     Reproduces the existing verification golds bit-for-bit.
 *   interface -- evaluate at the CO2-brine contact z_T + h, using the hydrostatic
 *     pressure p = pp_top + rho_n(pp_top)*|g|*h through the CO2 column. The
 *     thickness h = sat_n*H/(1-S_wr) is the sharp-interface form clamped to [0, H],
 *     with H = z_top - z_bottom from coupled geometry variables (NOT VEGeometry's
 *     ve_H), so interface mode works unchanged in FV models, which build H inline and
 *     do not run VEGeometry. h is built from sat_n here rather than from
 *     VEPlumeReconstruction's ve_h (which would form a dependency cycle, since
 *     VEPlumeReconstruction consumes ve_density). p carries AD wrt pp_top and sat_n.
 *
 * initQpStatefulProperties() seeds the old-value storage that
 * VEMassTimeDerivative / VEFVMassTimeDerivative read, so the storage term is correct
 * at t = 0. It uses the same reference depth as the regular compute: the geometry and
 * saturation come from coupled VARIABLE values, which are populated from the ICs
 * before stateful init, so the interface pressure is available there too.
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
  /// Fill ve_density/ve_viscosity at this qp. at_interface selects the EOS pressure:
  /// false -> pp_top (top surface); true -> the z_T + h interface pressure.
  void fillQpValues(bool at_interface);

  /// True if eos_reference_depth = interface (evaluate the EOS at z_T + h).
  const bool _eos_at_interface;

  /// Top-surface pore pressure [Pa] -- the top-surface EOS pressure (AD).
  const ADVariableValue & _pp_top;
  /// Constant formation temperature [K], used only to evaluate the EOS.
  const Real _temperature;

  /// Depth-averaged CO2 saturation, for the sharp-interface h (interface mode only).
  const ADVariableValue * const _sat_n;
  /// Residual water saturation in the CO2 zone [-] (interface mode only).
  const Real _S_wr;
  /// |g| [m/s2], for the hydrostatic increment (interface mode only).
  const Real _gravity_magnitude;
  /// Top/bottom elevations [m] giving H = z_top - z_bottom (interface mode only).
  const VariableValue * const _z_top;
  const VariableValue * const _z_bottom;

  const SinglePhaseFluidProperties & _fp_nw;
  const SinglePhaseFluidProperties & _fp_w;

  ADMaterialProperty<std::vector<Real>> & _density;
  ADMaterialProperty<std::vector<Real>> & _viscosity;
};
