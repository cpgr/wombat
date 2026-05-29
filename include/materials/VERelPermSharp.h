#pragma once

#include "Material.h"

/**
 * VERelPermSharp
 *
 * Upscaled relative permeabilities for a sharp CO2-brine interface
 * (Nordbotten and Celia, Geological Storage of CO2, Chapter 3).
 *
 * The CO2 plume occupies the top h metres of the formation. Below the
 * interface is fully brine-saturated rock. Depth-integrating the pointwise
 * relperm over each sub-column gives closed-form expressions in sat_n:
 *
 *   eff_sat  = sat_n / (1 - S_wr)          -- normalised saturation
 *   kr_n_up  = krn_max * eff_sat                -- CO2 upscaled relperm
 *   kr_w_up  = krw_max * (1 - eff_sat)          -- brine upscaled relperm
 *
 * where S_wr is the residual water saturation in the CO2 zone. eff_sat is
 * clamped to [0, 1] to handle numerical overshoot at the saturation front.
 *
 * These analytical curves are the first-class verification target: the
 * Nordbotten-Celia sloped-aquifer plume-nose velocity test depends on them.
 * They are also used in the 1D Buckley-Leverett test (flat formation, S_wr=0).
 *
 * For field cases with tabulated upscaled curves from fine-scale upscaling,
 * a separate VERelPermTabulated material will replace this one.
 *
 * Reads ve_saturation from VESaturation (index 0 = CO2).
 * Writes ve_relperm (index 0 = CO2, index 1 = brine).
 */
class VERelPermSharp : public Material
{
public:
  static InputParameters validParams();
  VERelPermSharp(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const Real _S_wr;     ///< Residual water saturation in the CO2 zone [-].
  const Real _krn_max;  ///< Maximum CO2 relative permeability [-].
  const Real _krw_max;  ///< Maximum brine relative permeability [-].

  /// Depth-averaged saturations from VESaturation [0]=CO2, [1]=brine.
  const ADMaterialProperty<std::vector<Real>> & _saturation;

  ADMaterialProperty<std::vector<Real>> & _relperm;
};
