#pragma once

#include "SinglePhaseFluidProperties.h"

/**
 * ConstantFluidProperties
 *
 * Minimal SinglePhaseFluidProperties returning a constant density and viscosity
 * that are independent of pressure and temperature (zero derivatives). The
 * FluidProperties module ships no constant-property class, so this provides the
 * verification-case counterpart of a full EOS.
 *
 * Its purpose is to let a SINGLE properties material (VEFluidProperties) be used
 * everywhere: ConstantFluidProperties reproduces the old VEFluidPropertiesConst
 * numbers exactly for the analytical verification golds (Buckley-Leverett,
 * Nordbotten-Celia, ...), while CO2FluidProperties / BrineFluidProperties (or any
 * other SinglePhaseFluidProperties) drive pressure/temperature-dependent field
 * cases through the same material -- no kernel changes, no re-golding.
 *
 * Only rho_from_p_T and mu_from_p_T are implemented (the two quantities the VE
 * flux and storage terms need). The AD overloads are generated automatically by
 * the SinglePhaseFluidProperties propfunc machinery from the Real derivative
 * forms, and since both derivatives are zero the AD result carries zero
 * derivatives -- correct for constant properties. Every other property inherits
 * the base-class "not implemented" error; if a future consumer needs one, add it
 * here.
 */
class ConstantFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();
  ConstantFluidProperties(const InputParameters & parameters);

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

  virtual std::string fluidName() const override { return "constant"; }

protected:
  /// Constant density [kg/m3].
  const Real _density;
  /// Constant dynamic viscosity [Pa*s].
  const Real _viscosity;
};
