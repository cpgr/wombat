#pragma once

#include "Material.h"

/**
 * VEFluidPropertiesConst
 *
 * Provides constant (pressure- and temperature-independent) fluid
 * properties for the two-phase CO2-brine system. Intended for
 * verification cases (Buckley-Leverett, Nordbotten-Celia) where a
 * simple, analytically tractable EOS is required.
 *
 * Declares:
 *   ve_density   [kg/m3]  -- ADMaterialProperty, size 2: [rho_CO2, rho_brine]
 *   ve_viscosity [Pa*s]   -- ADMaterialProperty, size 2: [mu_CO2,  mu_brine]
 *
 * Both properties are AD so the kernel interface is uniform with the
 * future pressure-dependent EOS (VEFluidStateBrineCO2). For constant
 * properties the AD derivatives are zero, which is correct.
 *
 * initQpStatefulProperties() seeds the old-value storage required by
 * VEMassTimeDerivative so the storage term is correct at t=0.
 */
class VEFluidPropertiesConst : public Material
{
public:
  static InputParameters validParams();
  VEFluidPropertiesConst(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

private:
  void fillQpValues();

  const Real _rho_co2;
  const Real _rho_brine;
  const Real _mu_co2;
  const Real _mu_brine;

  ADMaterialProperty<std::vector<Real>> & _density;
  ADMaterialProperty<std::vector<Real>> & _viscosity;
};
