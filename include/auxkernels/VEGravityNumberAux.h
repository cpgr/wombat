#pragma once

#include "AuxKernel.h"

/**
 * VEGravityNumberAux
 *
 * Computes the gravity number (Gamma), a dimensionless VE-validity indicator
 * that quantifies the ratio of vertical gravity segregation to lateral
 * viscous flow.  Large Gamma means fast vertical segregation relative to
 * horizontal flow -- the VE assumption is strong.  Small Gamma (< ~10)
 * indicates that vertical equilibrium may not have been reached and the
 * depth-integrated model should be interpreted with caution.
 *
 *   Gamma = k_v * delta_rho * g * H^2
 *           / ( mu_n * phi_bar * Q * L )
 *
 * where:
 *   k_v       = vertical permeability [m^2]
 *   delta_rho = rho_brine - rho_CO2  [kg/m^3]
 *   g         = gravitational acceleration [m/s^2]
 *   H         = formation thickness [m]  (spatial, from ve_H)
 *   mu_n      = CO2 viscosity [Pa*s]
 *   phi_bar   = depth-averaged porosity [-]  (spatial, from ve_phi_bar)
 *   Q         = characteristic volumetric injection rate [m^3/s]
 *   L         = characteristic horizontal length [m]
 *
 * H and phi_bar vary spatially so Gamma is computed per quadrature point
 * and averaged to the element.  k_v, delta_rho, mu_n, Q, L, and g are
 * uniform scalar parameters.
 *
 * Reads:
 *   - ve_H       : material property from VEGeometry
 *   - ve_phi_bar : material property from VEPorosity
 *
 * Parameters: k_v, delta_rho, mu_n, Q, L, gravity
 */
class VEGravityNumberAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VEGravityNumberAux(const InputParameters & parameters);

protected:
  Real computeValue() override;

  const MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _phi_bar;

  const Real _k_v;
  const Real _delta_rho;
  const Real _mu_n;
  const Real _Q;
  const Real _L;
  const Real _gravity_magnitude;
};
