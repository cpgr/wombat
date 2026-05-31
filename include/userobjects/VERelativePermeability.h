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
   * AD convenience: returns kr with its sat_n-derivative seeded from the
   * Real value + derivative. Used by both the FE and FV adapters.
   */
  ADReal relativePermeabilityAD(const ADReal & sat_n, unsigned int phase) const;
};
