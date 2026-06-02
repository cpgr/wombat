#pragma once

#include "SinglePhaseFluidProperties.h"

class BrineFluidProperties;

/**
 * SimpleBrineFluidProperties
 *
 * Fixed-salinity brine exposed through the two-variable (pressure, temperature)
 * SinglePhaseFluidProperties interface that VEFluidProperties consumes.
 *
 * The FluidProperties-module BrineFluidProperties is a three-variable EOS
 * (pressure, temperature, salt mass fraction), so it cannot be handed to
 * VEFluidProperties directly -- the VE interface carries only pressure and
 * temperature. This class bridges the gap: it holds a constant salt mass
 * fraction, constructs a BrineFluidProperties internally (the same way
 * BrineFluidProperties itself builds its component water/NaCl UserObjects), and
 * forwards rho/mu evaluations to it at the fixed salinity. The derivatives wrt
 * the third variable (salinity) are discarded, so only d/dp and d/dT are
 * exposed -- exactly what a constant-salinity model needs, and consistent with
 * the AD seeding done by SinglePhaseFluidProperties' propfunc machinery.
 *
 * Variable salinity is intentionally out of scope; if it is ever needed, brine
 * salinity would have to become a solution/coupled field and the VE interface
 * widened accordingly.
 *
 * Use as the wetting-phase EOS (fp_w) of VEFluidProperties when the wetting
 * phase is brine rather than pure water. Temperature is in Kelvin (as supplied
 * by VEFluidProperties' temperature parameter).
 *
 * Only rho_from_p_T and mu_from_p_T are implemented (what the VE flux/storage
 * terms need); other properties inherit the base-class "not implemented" error.
 */
class SimpleBrineFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();
  SimpleBrineFluidProperties(const InputParameters & parameters);

  // The Real overrides below would otherwise hide the base-class AD overloads;
  // bring them back into scope (they dispatch to the derivative forms here).
  using SinglePhaseFluidProperties::rho_from_p_T;
  using SinglePhaseFluidProperties::mu_from_p_T;

  virtual Real rho_from_p_T(Real p, Real T) const override;
  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;

  virtual Real mu_from_p_T(Real p, Real T) const override;
  virtual void
  mu_from_p_T(Real p, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const override;

  virtual std::string fluidName() const override { return "brine"; }

protected:
  /// Constant salt (NaCl) mass fraction [kg/kg].
  const Real _xnacl;
  /// Underlying three-variable brine EOS, constructed internally.
  const BrineFluidProperties * _brine_fp;
};
