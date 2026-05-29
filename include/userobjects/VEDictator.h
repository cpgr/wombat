#pragma once

#include "PorousFlowDictator.h"

/**
 * VEDictator -- the single source of truth for a VE simulation.
 *
 * Extends PorousFlowDictator with the VE flavour flag. Always configures
 * two fluid phases (0 = CO2/non-wetting, 1 = brine/wetting) and two
 * fluid components (0 = CO2, 1 = water). All VE Materials and Kernels
 * request a reference to this object so they can query configuration
 * consistently without duplicating parameter declarations.
 */
class VEDictator : public PorousFlowDictator
{
public:
  /// VE physics variant, controls which constitutive relations are active.
  enum class VEFlavour
  {
    SHARP_INTERFACE,  ///< Classical sharp CO2-brine interface, no capillary fringe.
    CAPILLARY_FRINGE, ///< Nordbotten-Dahle 2011 capillary fringe via lookup table.
    WITH_HYSTERESIS   ///< Capillary fringe + trapping via scanning curves.
  };

  static InputParameters validParams();
  VEDictator(const InputParameters & parameters);

  /// VE physics variant selected in the input file.
  VEFlavour flavour() const { return _flavour; }

  /// Phase index for the non-wetting (CO2) phase -- always 0.
  static constexpr unsigned int co2Phase() { return 0; }

  /// Phase index for the wetting (brine) phase -- always 1.
  static constexpr unsigned int brinePhase() { return 1; }

  /// Component index for CO2 -- always 0.
  static constexpr unsigned int co2Component() { return 0; }

  /// Component index for water -- always 1.
  static constexpr unsigned int waterComponent() { return 1; }

private:
  const VEFlavour _flavour;
};
