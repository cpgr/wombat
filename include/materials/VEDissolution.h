#pragma once

#include "Material.h"

/**
 * VEDissolution
 *
 * Convective-dissolution rate for CO2 dissolving into the underlying brine, as an areal
 * mass flux (per unit map-view area) consumed by VEDissolutionSink and accumulated by
 * VEDissolvedCO2Aux.
 *
 *   ve_dissolution_rate = q0 * gate(sat_n) * capacity(c_diss)   [kg / m^2 / s]
 *
 * q0 is the constant-flux convective dissolution rate. In the density-driven convective
 * regime (Pau et al. 2010; Neufeld et al. 2010; Hesse) the dissolution flux per unit
 * plume-footprint area becomes approximately constant (independent of H), q0 ~ alpha *
 * k * delta_rho_c * g / mu * C_sat. For v1 q0 is supplied directly as a parameter
 * (dissolution_flux); deriving it from the correlation is a later refinement.
 *
 * gate(sat_n) = min(1, sat_n / s_ref) rate-limits dissolution as the plume thins so the
 * sink cannot drive sat_n negative: full rate where the plume is present (sat_n >= s_ref),
 * tapering linearly to zero as sat_n -> 0 (no plume -> no dissolution).
 *
 * capacity = max(0, 1 - c_diss / c_cap) optionally stops dissolution once the column has
 * absorbed its CO2 capacity c_cap [kg/m^2] (e.g. rho_brine*phi*H*X_sat). c_diss is the
 * lagged areal dissolved mass (from VEDissolvedCO2Aux), so this factor is frozen within a
 * solve. Omit dissolved_co2 / c_cap to disable the cap (unbounded dissolution).
 *
 * Sink-only / Option A: the dissolved CO2 is tracked as an immobile areal inventory; it is
 * not added to the brine phase or transported. The freed pore volume is supplied by brine
 * inflow through the coupled pressure solve (open system).
 */
class VEDissolution : public Material
{
public:
  static InputParameters validParams();
  VEDissolution(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Constant-flux convective dissolution rate q0 [kg/m^2/s].
  const Real _q0;
  /// Gate reference saturation [-]: full rate for sat_n >= s_ref.
  const Real _s_ref;

  /// Whether the column-capacity cap is active.
  const bool _has_cap;
  /// Column CO2 capacity c_cap [kg/m^2]; only valid if _has_cap.
  const Real _c_cap;
  /// Lagged areal dissolved CO2 mass c_diss [kg/m^2] (from VEDissolvedCO2Aux).
  const VariableValue & _c_diss;

  /// Depth-averaged phase saturations [-] (AD); sat_n = [0].
  const ADMaterialProperty<std::vector<Real>> & _saturation;

  /// Output: areal dissolution rate [kg/m^2/s] (AD in sat_n via the gate).
  ADMaterialProperty<Real> & _rate;
};
