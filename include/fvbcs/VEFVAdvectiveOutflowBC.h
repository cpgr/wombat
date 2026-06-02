#pragma once

#include "FVFluxBC.h"
#include "VEDictator.h"

/**
 * VEFVAdvectiveOutflowBC
 *
 * Open / advective-outflow boundary condition for the depth-integrated VE flow,
 * the boundary counterpart of VEFVAdvectiveFlux. It lets a phase leave the domain
 * at the Darcy rate instead of piling up against a no-flow wall (the default
 * natural BC). Assign it to variable = sat_n (fluid_phase = 0) on the updip /
 * downstream boundary so migrating CO2 can exit; the brine equation is usually
 * handled by a Dirichlet pp_top there.
 *
 *   F_c . n = -H * K_nn * (kr_c * rho_c / mu_c)
 *                       * ( grad(pp_top).n + rho_c * |g| * grad(z_T).n )
 *
 * evaluated on the boundary face with the SAME functors as the interior kernel.
 * The pressure drive picks up the pp_top Dirichlet BC (so the boundary Darcy
 * velocity is set by the prescribed pressure plus buoyancy), and the relative
 * permeability is taken from the interior cell -- the upwind side for outflow --
 * so the flux carries the mobility of the phase actually present inside.
 *
 * It is strictly OUTFLOW: when the boundary Darcy velocity points inward (or is
 * zero) the flux is set to zero, since the exterior holds none of this phase.
 * Thus it never draws the phase back into the domain -- it only lets a phase that
 * reaches the boundary leave at the Darcy rate. (Brine outflow/inflow is normally
 * handled by a Dirichlet pp_top on the same boundary.)
 *
 * Capillary and TVD options are intentionally omitted: an open far boundary does
 * not need front sharpening, and grad(Pc^up) there is negligible.
 */
class VEFVAdvectiveOutflowBC : public FVFluxBC
{
public:
  static InputParameters validParams();
  VEFVAdvectiveOutflowBC(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const VEDictator & _dictator;
  const unsigned int _fluid_phase;
  const Real _gravity_magnitude;

  // --- Functors (boundary-aware) ---
  const Moose::Functor<ADReal> & _pp_top;
  const Moose::Functor<ADReal> & _z_top;
  const Moose::Functor<ADReal> & _z_bottom;
  const Moose::Functor<ADReal> & _relperm; ///< ve_relperm_n (phase 0) or ve_relperm_w (phase 1)

  // --- Elem-side material properties (constant; valid on boundary faces) ---
  const MaterialProperty<RealTensorValue> & _K_up;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _viscosity;
};
