#pragma once

#include "GeneralUserObject.h"
#include "MooseTypes.h"

/**
 * VERelativePermeability
 *
 * Abstract base class for VE upscaled relative-permeability models, exposed as a
 * swappable UserObject (the same pattern as PorousFlowCapillaryPressure /
 * VECapillaryPressureTable for capillary pressure).
 *
 * A model implements the relative permeability and its derivative as functions
 * of the depth-averaged CO2 saturation, returning plain Real (value + derivative,
 * following the PorousFlowCapillaryPressure convention):
 *
 *   relativePermeability (sat_n, phase)
 *   dRelativePermeability(sat_n, phase)   = d(kr)/d(sat_n)
 *
 * with phase 0 = CO2 (non-wetting) and phase 1 = brine (wetting).
 *
 * Both consumers are thin adapters that delegate the physics here:
 *   VERelPerm    (Material)        -> ve_relperm   qp property, for FE kernels
 *   VEFVRelPerm  (FunctorMaterial) -> ve_relperm_n/w functors, for the FV kernel
 *
 * The non-virtual relativePermeabilityAD() assembles the AD value (seeding the
 * sat_n derivative from value + derivative) so neither adapter has to repeat that.
 *
 * Hysteresis (residual trapping) is supported through history-aware overloads that
 * additionally take sat_n_max, the maximum historical CO2 saturation (the turning
 * point of the scanning curve, carried by a lagged stateful aux -- see
 * VESaturationMaxAux). The base-class defaults IGNORE sat_n_max and delegate to the
 * two-argument (drainage) form, so a non-hysteretic model is automatically a no-op on
 * history and every existing subclass / input is unchanged. A hysteretic model
 * (VERelPermHysteresisUO) overrides the three-argument forms and trappedSaturation().
 * sat_n_max is FROZEN within a solve (lagged aux), so it enters as a plain Real and
 * the AD seed remains purely d(kr)/d(sat_n).
 */
class VERelativePermeability : public GeneralUserObject
{
public:
  static InputParameters validParams();
  VERelativePermeability(const InputParameters & parameters);

  void initialize() final {}
  void execute() final {}
  void finalize() final {}

  /// Relative permeability of @p phase (0 = CO2, 1 = brine) at CO2 saturation @p sat_n.
  virtual Real relativePermeability(Real sat_n, unsigned int phase) const = 0;

  /// Derivative d(relativePermeability)/d(sat_n).
  virtual Real dRelativePermeability(Real sat_n, unsigned int phase) const = 0;

  /**
   * History-aware relative permeability on the scanning curve through turning point
   * @p sat_n_max. Default ignores history (drainage form); overridden by hysteretic
   * models.
   */
  virtual Real relativePermeability(Real sat_n, Real /*sat_n_max*/, unsigned int phase) const
  {
    return relativePermeability(sat_n, phase);
  }

  /// History-aware d(kr)/d(sat_n). Default ignores history (drainage form).
  virtual Real dRelativePermeability(Real sat_n, Real /*sat_n_max*/, unsigned int phase) const
  {
    return dRelativePermeability(sat_n, phase);
  }

  /**
   * Depth-averaged residually trapped CO2 saturation given the current saturation and
   * the turning point @p sat_n_max. Default is zero (no trapping); hysteretic models
   * return the Land/Killough trapped saturation. Consumed by VETrappedSaturationAux
   * (-> the trapped/mobile CO2 mass postprocessors).
   */
  virtual Real trappedSaturation(Real /*sat_n*/, Real /*sat_n_max*/) const { return 0.0; }

  /**
   * AD convenience: returns kr with its sat_n-derivative seeded from the
   * Real value + derivative. Used by both the FE and FV adapters.
   */
  ADReal relativePermeabilityAD(const ADReal & sat_n, unsigned int phase) const;

  /// History-aware AD overload (sat_n_max frozen): seeds the derivative wrt sat_n only.
  ADReal
  relativePermeabilityAD(const ADReal & sat_n, Real sat_n_max, unsigned int phase) const;
};
