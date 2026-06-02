#pragma once

#include "FunctorMaterial.h"
#include "SinglePhaseFluidProperties.h"

/**
 * VEFVFluidProperties
 *
 * FV adapter: publishes the phase density and viscosity as the functors
 * ve_density_n / ve_density_w / ve_viscosity_n / ve_viscosity_w for the FV flux
 * kernel, delegating to two SinglePhaseFluidProperties UserObjects (fp_nw / fp_w)
 * evaluated at a reference pressure derived from the pp_top functor and a constant
 * temperature.
 *
 * The functors evaluate pp_top at whatever face argument the flux kernel supplies.
 * VEFVAdvectiveFlux evaluates them at a central-difference (face-averaged) argument,
 * so the density and viscosity used in the depth-integrated Darcy flux are
 * face-correct -- they depend on the pressure of BOTH adjacent cells and carry AD
 * wrt pp_top on both sides of the face. This is the FV counterpart of the elemental
 * VEFluidProperties material; both wrap the same FP UserObjects at the same
 * temperature and the same reference-depth strategy, so they are consistent by
 * construction.
 *
 * EOS reference depth (CLAUDE.md key-subtlety 2), switchable via
 * eos_reference_depth, mirrors VEFluidProperties:
 *   top_surface (default) -- evaluate at the pp_top functor (the formation top z_T).
 *   interface -- evaluate at the CO2-brine contact z_T + h, using the hydrostatic
 *     pressure p = pp_top + rho_n(pp_top)*|g|*h with the sharp-interface thickness
 *     h = sat_n*H/(1-S_wr) (H = z_top - z_bottom) clamped to [0, H], all evaluated at
 *     the same face argument so the reference pressure stays face-correct.
 *
 * VEFluidProperties is still required in an FV input: the FV mass-storage kernel
 * (VEFVMassTimeDerivative) is elemental and reads the stateful ve_density (with its
 * old value) from VEFluidProperties. This functor adapter feeds only the flux kernel.
 */
class VEFVFluidProperties : public FunctorMaterial
{
public:
  static InputParameters validParams();
  VEFVFluidProperties(const InputParameters & parameters);

protected:
  /// EOS reference pressure at the given functor argument: pp_top in top_surface mode,
  /// or the z_T + h interface pressure in interface mode. Templated on the functor
  /// spatial/state argument types so it works for any face/elem evaluation.
  template <typename Space, typename State>
  ADReal refPressure(const Space & r, const State & t) const
  {
    const ADReal p_top = _pp_top(r, t);
    if (!_eos_at_interface)
      return p_top;

    const ADReal H = (*_z_top)(r, t) - (*_z_bottom)(r, t);
    ADReal h = (*_sat_n)(r, t) * H / (1.0 - _S_wr);
    if (MetaPhysicL::raw_value(h) <= 0.0)
      h = 0.0;
    else if (MetaPhysicL::raw_value(h) >= MetaPhysicL::raw_value(H))
      h = H;

    const ADReal rho_n_top = _fp_nw.rho_from_p_T(p_top, ADReal(_temperature));
    return p_top + rho_n_top * _gravity_magnitude * h;
  }

  /// True if eos_reference_depth = interface (evaluate the EOS at z_T + h).
  const bool _eos_at_interface;

  /// Top-surface pore pressure functor [Pa] -- the top-surface EOS pressure.
  const Moose::Functor<ADReal> & _pp_top;
  /// Constant formation temperature [K], used only to evaluate the EOS.
  const Real _temperature;

  /// Residual water saturation in the CO2 zone [-] (used in interface mode).
  const Real _S_wr;
  /// |g| [m/s2], for the hydrostatic increment (used in interface mode).
  const Real _gravity_magnitude;
  /// Functors for the sharp-interface h = sat_n*(z_top - z_bottom)/(1 - S_wr); only
  /// bound in interface mode (nullptr otherwise).
  const Moose::Functor<ADReal> * const _sat_n;
  const Moose::Functor<ADReal> * const _z_top;
  const Moose::Functor<ADReal> * const _z_bottom;

  const SinglePhaseFluidProperties & _fp_nw;
  const SinglePhaseFluidProperties & _fp_w;
};
