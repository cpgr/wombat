#pragma once

#include "VERelativePermeability.h"

/**
 * VERelPermHysteresisUO
 *
 * Hysteretic upscaled relative permeability with residual (capillary) trapping,
 * the saturation-space sharp-interface analogue of the Land/Killough model. The
 * non-wetting (CO2) curve depends on the turning point sat_n_max (max historical
 * CO2 saturation, supplied by VESaturationMaxAux); the wetting (brine) curve is
 * taken non-hysteretic (a common, defensible simplification).
 *
 * In effective CO2 saturation s = sat_n / (1 - S_wr), s_max = sat_n_max / (1 - S_wr),
 * with Land coefficient C such that the maximum trapped CO2 saturation after full
 * drainage (s_max = 1) is S_gr_max:
 *
 *   C = (1 - S_wr)/S_gr_max - 1      (so s_gr_max = S_gr_max/(1-S_wr) = 1/(1+C))
 *
 * Drainage (sat_n >= sat_n_max):    kr_n = krn_max * s            (all CO2 mobile)
 * Imbibition (sat_n < sat_n_max):   s_gr = s_max / (1 + C*s_max)
 *                                   kr_n = krn_max * s_max * (s - s_gr)/(s_max - s_gr)
 *
 * The imbibition line is continuous with drainage at the turning point
 * (kr_n = krn_max*s_max at s = s_max) and vanishes at s = s_gr, immobilising the
 * residual CO2. The trapped saturation is defined consistently with this kr so that
 * the mobile saturation sat_n - sat_n_trap equals the kr-implied mobile saturation:
 *
 *   sat_n_trap = (1 - S_wr) * s_gr * (s_max - s)/(s_max - s_gr)
 *
 * which is 0 at the turning point (continuous reversal), grows monotonically as the
 * column imbibes, and reaches the Land residual (1 - S_wr)*s_gr when s = s_gr.
 *
 * Pure drainage (monotonically increasing sat_n) never enters the imbibition branch,
 * so this UO reduces EXACTLY to VERelPermSharpUO during injection -- a hysteresis run
 * with monotonic loading reproduces the sharp-interface result.
 *
 * kr_w (brine) = krw_max * (1 - s) on both branches: trapped CO2 still occupies pore
 * space and blocks brine, so the wetting mobility depends on the total CO2 saturation.
 */
class VERelPermHysteresisUO : public VERelativePermeability
{
public:
  static InputParameters validParams();
  VERelPermHysteresisUO(const InputParameters & parameters);

  // Two-argument (drainage / no-history) forms.
  virtual Real relativePermeability(Real sat_n, unsigned int phase) const override;
  virtual Real dRelativePermeability(Real sat_n, unsigned int phase) const override;

  // History-aware (scanning-curve) forms.
  virtual Real
  relativePermeability(Real sat_n, Real sat_n_max, unsigned int phase) const override;
  virtual Real
  dRelativePermeability(Real sat_n, Real sat_n_max, unsigned int phase) const override;
  virtual Real trappedSaturation(Real sat_n, Real sat_n_max) const override;

protected:
  /// Effective CO2 saturation s = clamp(sat_n/(1-S_wr), 0, 1).
  Real effSat(Real sat_n) const;

  const Real _S_wr;
  const Real _krn_max;
  const Real _krw_max;
  const Real _S_gr_max; ///< Max trapped CO2 saturation after full drainage [-].
  const Real _C;        ///< Land coefficient, C = (1-S_wr)/S_gr_max - 1.
};
