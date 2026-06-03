#pragma once

#include "FunctorMaterial.h"

/**
 * VEFVCapPressure
 *
 * FV functor material: declares the functor ve_pc_up for use by VEFVAdvectiveFlux
 * when capillary=true. The functor is evaluated at any face argument, so the FV
 * kernel can call .gradient(face, state) to obtain the face-normal gradient of
 * the upscaled capillary pressure via the standard MOOSE FV reconstruction.
 *
 * Sharp-interface formula (Nordbotten-Celia), with delta_rho = rho_w - rho_n read
 * per evaluation from the density functors (so it tracks pressure/temperature):
 *   Pc^up(r,t)        = delta_rho * |g| * sat_n(r,t) * H(r,t) / (1 - S_wr) + pc_entry
 *   ve_dPcup_dsatn    = delta_rho * |g| * H(r,t) / (1 - S_wr)
 *
 * where H(r,t) = z_top(r,t) - z_bottom(r,t).
 *
 * VEFVAdvectiveFlux consumes ve_dPcup_dsatn (not ve_pc_up) and forms
 * grad(Pc^up).n = ve_dPcup_dsatn(face) * grad(sat_n).n, using the sat_n VARIABLE
 * gradient (boundary-aware). ve_pc_up is published for diagnostics only: taking
 * the Green-Gauss gradient of a material functor does not pick up the sat_n
 * Dirichlet BC at a boundary face, so it must not be used for the flux.
 *
 * Parameters:
 *   sat_n      -- name of the FV saturation functor (typically the variable name)
 *   z_top      -- name of the z_T functor (FV aux variable or variable)
 *   z_bottom   -- name of the z_B functor
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
  const Real _S_wr;
  const Real _gravity_magnitude;
  const Real _pc_entry;

  const Moose::Functor<ADReal> & _sat_n;
  const Moose::Functor<ADReal> & _z_top;
  const Moose::Functor<ADReal> & _z_bottom;
  const Moose::Functor<ADReal> & _rho_n;
  const Moose::Functor<ADReal> & _rho_w;
};
