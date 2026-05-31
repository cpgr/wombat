#pragma once

#include "Material.h"

/**
 * VEUpscaledCapPressure
 *
 * Computes the upscaled (macroscopic) capillary pressure of the depth-integrated
 * column as a function of the reconstructed plume thickness:
 *
 *   Pc^up = (rho_w - rho_n) * |g| * h  +  pc_entry
 *
 * The (rho_w - rho_n) * |g| * h term is the buoyant pressure difference across a
 * CO2 column of height h (sharp-interface upscaled Pc; Nordbotten & Celia).
 *
 * Also declares ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = (rho_w-rho_n)*|g|*H/(1-S_wr)
 * (sharp-interface formula), consumed by VEAdvectiveFluxS when capillary=true to
 * form grad(Pc^up) = ve_dPcup_dsatn * grad(sat_n) at each quadrature point.
 *
 * Both ve_pc_up and ve_dPcup_dsatn are AD properties so the Jacobian propagates
 * through density (when density becomes pressure-dependent in a future upgrade).
 *
 * Reads:
 *   - ve_h       : ADMaterialProperty<Real> from VEPlumeReconstruction
 *   - ve_density : ADMaterialProperty<std::vector<Real>> from VEFluidProperties
 *   - ve_H       : MaterialProperty<Real> from VEGeometry (for ve_dPcup_dsatn)
 */
class VEUpscaledCapPressure : public Material
{
public:
  static InputParameters validParams();
  VEUpscaledCapPressure(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const Real _gravity_magnitude;
  const Real _pc_entry;
  const Real _S_wr;

  const ADMaterialProperty<Real> & _h;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const MaterialProperty<Real> & _H;

  ADMaterialProperty<Real> & _pc_up;
  MaterialProperty<Real> & _dPcup_dsatn;
};
