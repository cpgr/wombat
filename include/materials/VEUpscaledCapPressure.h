#pragma once

#include "Material.h"

class PorousFlowCapillaryPressure;

/**
 * VEUpscaledCapPressure
 *
 * Computes the upscaled (macroscopic) capillary pressure of the depth-integrated
 * column as a function of the reconstructed plume thickness:
 *
 *   Pc^up = (rho_w - rho_n) * |g| * h  +  pc_entry
 *
 * The (rho_w - rho_n) * |g| * h term is the buoyant pressure difference across a
 * CO2 column of height h (Nordbotten & Celia / Nordbotten & Dahle).
 *
 * Also declares the two partial derivatives of Pc^up consumed by VEAdvectiveFluxS
 * when capillary != none, to assemble the full chain-rule gradient
 *
 *   grad(Pc^up) = ve_dPcup_dsatn * grad(sat_n) + ve_dPcup_dH * grad(H)
 *               = (rho_w - rho_n) * |g| * grad(h)
 *
 * at each quadrature point. Because grad(Pc^up) = delta_rho*g*grad(h) in BOTH
 * reconstruction modes, the coefficients are just the chain-rule partials of
 * h(sat_n, H):
 *
 *   ve_dPcup_dsatn = d(Pc^up)/d(sat_n) = delta_rho*|g| * dh/dsat_n
 *   ve_dPcup_dH    = d(Pc^up)/d(H)     = delta_rho*|g| * dh/dH
 *
 * with, writing S_n(h) = 1 - Sw(delta_rho*|g|*h) for the CO2 saturation at the
 * plume top,
 *
 *   sharp_interface:  dh/dsat_n = H/(1-S_wr),  dh/dH = sat_n/(1-S_wr)
 *   capillary_fringe: dh/dsat_n = H/S_n(h),    dh/dH = sat_n/S_n(h)
 *
 * Sharp interface is the special case S_n(h) = 1 - S_wr (constant).
 *
 * ve_dPcup_dsatn is an AD property: in sharp_interface mode it is constant in
 * sat_n (zero AD derivative for constant density, so the kernel's grad(sat_n)
 * term is unchanged from the older plain-Real coefficient); in capillary_fringe
 * mode it depends on sat_n through S_n(h), and the AD derivative -- propagated
 * from the AD ve_h via the PorousFlowCapillaryPressure AD saturation overload --
 * makes the kernel's grad(sat_n) Jacobian exact. ve_dPcup_dH is AD because
 * grad(H) is fixed geometry (carries no Jacobian), so the d/d(sat_n) derivative
 * of the grad(H) term must ride on this coefficient.
 *
 * Reads:
 *   - ve_h       : ADMaterialProperty<Real> from VEPlumeReconstruction
 *   - ve_density : ADMaterialProperty<std::vector<Real>> from VEFluidProperties
 *   - ve_H       : MaterialProperty<Real> from VEGeometry
 *   - sat_n      : coupled primary variable (capillary_fringe mode only)
 */
class VEUpscaledCapPressure : public Material
{
public:
  static InputParameters validParams();
  VEUpscaledCapPressure(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const MooseEnum _mode;
  const Real _gravity_magnitude;
  const Real _pc_entry;
  const Real _S_wr;

  /// CO2 saturation at the plume top, S_n(h) = 1 - Sw(delta_rho*|g|*h), used in
  /// the capillary_fringe partials. Floors near-zero values (plume nose) to keep
  /// the H/S_n divisions finite, mirroring VEPlumeReconstruction.
  ADReal saturationAtPlumeTop(const ADReal & delta_rho) const;

  const ADMaterialProperty<Real> & _h;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const MaterialProperty<Real> & _H;

  // capillary_fringe only
  const PorousFlowCapillaryPressure * _pc_uo;
  const ADVariableValue & _sat_n;

  ADMaterialProperty<Real> & _pc_up;
  ADMaterialProperty<Real> & _dPcup_dsatn;
  ADMaterialProperty<Real> & _dPcup_dH;
};
