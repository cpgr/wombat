#pragma once

#include "Material.h"

/**
 * VESaturation
 *
 * Wraps the primary variable sat_n (depth-averaged CO2 saturation)
 * into the two-element saturation vector expected by the VE kernels:
 *
 *   ve_saturation[0] = sat_n        (CO2, non-wetting)
 *   ve_saturation[1] = 1 - sat_n    (brine, wetting)
 *
 * The property is declared as ADMaterialProperty so that AD derivatives
 * of sat_n propagate into VEMassTimeDerivative and VERelPerm
 * without any hand-coded chain rule.
 *
 * initQpStatefulProperties() seeds the old-value storage from the initial
 * condition so that VEMassTimeDerivative computes a correct storage term
 * at the first timestep.
 */
class VESaturation : public Material
{
public:
  static InputParameters validParams();
  VESaturation(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// AD value of sat_n -- carries Jacobian derivatives through to consumers.
  const ADVariableValue & _sat_n;

  /// Non-AD value used only in initQpStatefulProperties for the initial seed.
  const VariableValue & _sat_n_real;

  ADMaterialProperty<std::vector<Real>> & _saturation;
};
