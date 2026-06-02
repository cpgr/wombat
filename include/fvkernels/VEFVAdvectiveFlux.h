#pragma once

#include "FVFluxKernel.h"
#include "VEDictator.h"

/**
 * VEFVAdvectiveFlux
 *
 * Finite-volume depth-integrated Darcy advective flux for one fluid phase:
 *
 *   F_c . n = -H * K_nn * (kr_c * rho_c / mu_c)
 *                       * ( grad(pp_top).n + rho_c * |g| * grad(z_T).n )
 *
 * A single kernel serves both equations -- assign it to variable = sat_n with
 * fluid_phase = 0 (CO2) and to variable = pp_top with fluid_phase = 1 (brine).
 *
 * Everything that varies across a face is evaluated through the MOOSE functor
 * system at a single-sided face argument, which is boundary-aware (uses the
 * face's BC value on a boundary face, reconstruction internally):
 *
 *   H           = z_top(face) - z_bottom(face)
 *   grad(z_T).n = z_top.gradient(face) . n
 *   grad(pp).n  = pp_top.gradient(face) . n            (carries AD wrt pp_top)
 *   kr_c        = ve_relperm_{n,w}( makeFace(upwind) )  (carries AD wrt sat_n)
 *   rho_c, mu_c = ve_{density,viscosity}_{n,w}(face)   (face avg; AD wrt pp_top)
 *
 * The relative permeability is the advected (upwinded) quantity: it is evaluated
 * at the UPWIND face via makeFace, so on an inflow boundary it picks up the
 * Dirichlet sat_n value (driving inflow correctly) and on internal faces it
 * selects the upstream cell -- the standard FVMatAdvection idiom. Because all
 * face quantities come from functors, the kernel never indexes a neighbor
 * material/array and is therefore safe on boundary faces (no neighbor element).
 *
 * rho_c and mu_c are evaluated through functors (ve_density_{n,w} /
 * ve_viscosity_{n,w} from VEFVFluidProperties) at a CENTRAL-DIFFERENCE
 * (face-averaged) argument, so for a pressure/temperature-dependent EOS they are
 * face-correct: they use both adjacent cell pressures and carry AD wrt pp_top on
 * both sides. A face average (rather than upwinding rho_c) keeps the smooth
 * density out of the monotone-limited advection and avoids the circular dependence
 * between the upwind direction and the buoyancy term. For constant fluid properties
 * the face average reduces to that constant, so this is exact and changes nothing.
 * (K is still taken elem-side; revisit harmonic K averaging when heterogeneous K is
 * introduced.)
 */
class VEFVAdvectiveFlux : public FVFluxKernel
{
public:
  static InputParameters validParams();
  VEFVAdvectiveFlux(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const VEDictator & _dictator;
  const unsigned int _fluid_phase;
  const Real _gravity_magnitude;

  /// If true, adds grad(Pc^up).n = ve_dPcup_dsatn*grad(sat_n).n +
  /// ve_dPcup_dH*grad(H).n to the CO2 (phase 0) Darcy potential. Requires
  /// VEFVCapPressure to be present (declares ve_dPcup_dsatn, ve_dPcup_dH).
  const bool _capillary;

  /// Face interpolation scheme for the advected relperm (default Upwind).
  Moose::FV::InterpMethod _advected_interp_method;

  // --- Functors (boundary-aware) ---
  const Moose::Functor<ADReal> & _pp_top;
  const Moose::Functor<ADReal> & _z_top;
  const Moose::Functor<ADReal> & _z_bottom;
  const Moose::Functor<ADReal> & _relperm; ///< ve_relperm_n (phase 0) or ve_relperm_w (phase 1)
  /// ve_density_n (phase 0) or ve_density_w (phase 1), evaluated face-averaged.
  const Moose::Functor<ADReal> & _density;
  /// ve_viscosity_n (phase 0) or ve_viscosity_w (phase 1), evaluated face-averaged.
  const Moose::Functor<ADReal> & _viscosity;
  /// d(Pc^up)/d(sat_n) coefficient functor (from VEFVCapPressure); null when
  /// capillary=false. Multiplied by grad(sat_n).n to form grad(Pc^up).n.
  const Moose::Functor<ADReal> * const _dPcup_dsatn;
  /// d(Pc^up)/d(H) coefficient functor (from VEFVCapPressure); null when
  /// capillary=false. Multiplied by grad(H).n to add the topography part of
  /// grad(Pc^up).n.
  const Moose::Functor<ADReal> * const _dPcup_dH;

  // --- Elem-side material property (constant; valid on boundary faces) ---
  const MaterialProperty<RealTensorValue> & _K_up;
};
