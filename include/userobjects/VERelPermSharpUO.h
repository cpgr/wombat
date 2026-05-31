#pragma once

#include "VERelativePermeability.h"

/**
 * VERelPermSharpUO
 *
 * Sharp-interface (Nordbotten-Celia) upscaled relative permeability:
 *
 *   s_eff = clamp( sat_n / (1 - S_wr), 0, 1 )
 *   kr_n  = krn_max * s_eff           (phase 0, CO2)
 *   kr_w  = krw_max * (1 - s_eff)     (phase 1, brine)
 *
 * This is the analytical verification curve. Tabulated / hysteretic curves from
 * the upscaling workflow are implemented as separate VERelativePermeability
 * subclasses; the FE (VERelPerm) and FV (VEFVRelPerm) adapters are unchanged.
 */
class VERelPermSharpUO : public VERelativePermeability
{
public:
  static InputParameters validParams();
  VERelPermSharpUO(const InputParameters & parameters);

  Real relativePermeability(Real sat_n, unsigned int phase) const override;
  Real dRelativePermeability(Real sat_n, unsigned int phase) const override;

protected:
  const Real _S_wr;    ///< Residual water saturation in the CO2 zone [-].
  const Real _krn_max; ///< Maximum CO2 (non-wetting) relative permeability [-].
  const Real _krw_max; ///< Maximum brine (wetting) relative permeability [-].
};
