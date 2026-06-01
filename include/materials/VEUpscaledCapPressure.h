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
 * Also declares the two partial derivatives of Pc^up (sharp-interface formulae),
 * consumed by VEAdvectiveFluxS when capillary=true to assemble the full
 * chain-rule gradient
 *
 *   grad(Pc^up) = ve_dPcup_dsatn * grad(sat_n) + ve_dPcup_dH * grad(H)
 *
 * at each quadrature point:
 *
 *   ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = (rho_w-rho_n)*|g|*H/(1-S_wr)
 *   ve_dPcup_dH    = d(Pc^up)/d(H)     = (rho_w-rho_n)*|g|*sat_n/(1-S_wr)
 *
 * ve_dPcup_dsatn is plain Real: the grad(sat_n) term's Jacobian rides on the
 * AD variable gradient _grad_u inside the kernel. ve_dPcup_dH is AD because
 * grad(H) is fixed geometry (carries no Jacobian), so the d/d(sat_n) derivative
 * of the grad(H) term must come from this coefficient. ve_pc_up is AD so the
 * Jacobian propagates through density (when density becomes pressure-dependent
 * in a future upgrade).
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
  ADMaterialProperty<Real> & _dPcup_dH;
};
