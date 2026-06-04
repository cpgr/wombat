#pragma once

#include "FunctorMaterial.h"

class PorousFlowCapillaryPressure;

/**
 * VEFVCapPressure
 *
 * FV functor material: declares the coefficient functors ve_dPcup_dsatn and
 * ve_dPcup_dH (and the diagnostic ve_pc_up) consumed by VEFVAdvectiveFlux when
 * capillary != none. The kernel forms grad(Pc^up).n = ve_dPcup_dsatn(face) *
 * grad(sat_n).n + ve_dPcup_dH(face) * grad(H).n, using the sat_n VARIABLE
 * gradient (boundary-aware) -- a material-functor Green-Gauss gradient would
 * miss the sat_n Dirichlet BC, so ve_pc_up is published for diagnostics only.
 *
 * Both modes share grad(Pc^up) = delta_rho*|g|*grad(h), so the coefficients are
 * the chain-rule partials of h(sat_n, H), with delta_rho = rho_w - rho_n read
 * per evaluation from the density functors and S_n(h) = 1 - Sw(delta_rho*|g|*h):
 *
 *   sharp_interface:  h = sat_n*H/(1-S_wr)
 *                     ve_dPcup_dsatn = delta_rho*|g|*H/(1-S_wr)
 *                     ve_dPcup_dH    = delta_rho*|g|*sat_n/(1-S_wr)
 *   capillary_fringe: h = Newton-inverted column equilibrium (Nordbotten & Dahle)
 *                     ve_dPcup_dsatn = delta_rho*|g|*H/S_n(h)
 *                     ve_dPcup_dH    = delta_rho*|g|*sat_n/S_n(h)
 *
 * where H = z_top - z_bottom. The FV flux Jacobian is AD-consistent in both
 * modes: the coefficient functor value carries d/d(sat_n) (in fringe via the
 * inverse-function-theorem AD seed on h), and it multiplies the AD gradUDotNormal.
 *
 * Parameters:
 *   sat_n      -- name of the FV saturation functor (typically the variable name)
 *   z_top      -- name of the z_T functor (FV aux variable or variable)
 *   z_bottom   -- name of the z_B functor
 *   mode       -- sharp_interface (default) or capillary_fringe
 *   pc_uo      -- PorousFlowCapillaryPressure UO (required in capillary_fringe)
 *   S_wr       -- residual water saturation [-]; must match VEPlumeReconstruction
 *   density_nw -- non-wetting (CO2) density functor [kg/m3] (default ve_density_n)
 *   density_w  -- wetting (brine) density functor [kg/m3] (default ve_density_w)
 *   gravity    -- gravity vector [m/s2]; only the magnitude is used
 *   pc_entry   -- constant entry pressure offset [Pa] (default 0)
 *
 * FE counterpart: VEUpscaledCapPressure (qp-based material).
 */
class VEFVCapPressure : public FunctorMaterial
{
public:
  static InputParameters validParams();
  VEFVCapPressure(const InputParameters & parameters);

protected:
  /// Depth-averaged saturation for a plume thickness h, column height H, and
  /// density difference delta_rho (rho_brine - rho_CO2) -- the capillary_fringe
  /// forward map inverted by plumeThickness (trapezoidal, 64 intervals).
  Real satNBar(Real h, Real H, Real delta_rho) const;

  /// Plume thickness h(sat_n) as an AD value. sharp_interface: closed form.
  /// capillary_fringe: Newton inversion on raw values, AD seeded via the
  /// inverse-function-theorem dh/dsat_n = H/S_n(h) (mirrors VEPlumeReconstruction).
  ADReal plumeThickness(const ADReal & delta_rho, const ADReal & H, const ADReal & sat_n) const;

  /// S_n(h) = 1 - Sw(delta_rho*|g|*h) at the plume top, floored near zero.
  ADReal saturationTop(const ADReal & delta_rho, const ADReal & h) const;

  const MooseEnum _mode;
  const Real _S_wr;
  const Real _gravity_magnitude;
  const Real _pc_entry;

  const Moose::Functor<ADReal> & _sat_n;
  const Moose::Functor<ADReal> & _z_top;
  const Moose::Functor<ADReal> & _z_bottom;
  const Moose::Functor<ADReal> & _rho_n;
  const Moose::Functor<ADReal> & _rho_w;

  // capillary_fringe only
  const PorousFlowCapillaryPressure * _pc_uo;
};
