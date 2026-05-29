#pragma once

#include "Material.h"

class PorousFlowCapillaryPressure;

/**
 * VEPlumeReconstruction
 *
 * Reconstructs the local CO2 plume thickness ve_h [m] from the depth-averaged
 * saturation sat_n. Two modes are supported:
 *
 *   sharp_interface (default):
 *     Nordbotten-Celia analytical formula:
 *       h = sat_n * H / (1 - S_wr)
 *
 *   capillary_fringe (Nordbotten & Dahle 2011):
 *     Newton inversion of the depth-averaged saturation integral:
 *       sat_n_bar(h) = (1/H) * integral_0^h [1 - Sw(Pc(zeta))] dzeta
 *     where Pc(zeta) = (rho_brine - rho_CO2) * g * zeta is the buoyancy-driven
 *     capillary pressure at height zeta above the CO2-brine contact. The density
 *     difference is read per quadrature point from ve_density (supplied by a
 *     VEFluidProperties material), and Sw(Pc) is supplied by a
 *     PorousFlowCapillaryPressure UserObject.
 *     The Newton derivative dF/dh = S_n(h)/H follows from the Leibniz rule.
 *     The integral is evaluated by the trapezoidal rule (64 intervals).
 *     The sharp-interface estimate is used as the warm start.
 *
 * S_wr is required in both modes: used directly in the sharp_interface formula
 * and implicit in the capillary_fringe Pc curve (declared here for consistency
 * with downstream materials such as VERelPermSharp).
 *
 * The result is clamped to [0, H] in both modes to absorb numerical overshoot.
 *
 * Declares the material property ve_h, consumed by VEUpscaledRelPerm,
 * VEUpscaledCapPressure, and VEFluidState. VEPlumeHeightAux reads ve_h for output.
 *
 * Reads:
 *   - sat_n : coupled primary variable (depth-averaged CO2 saturation)
 *   - ve_H  : material property from VEGeometry (formation thickness [m])
 */
class VEPlumeReconstruction : public Material
{
public:
  static InputParameters validParams();
  VEPlumeReconstruction(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Depth-averaged saturation for a given plume thickness h, column height H,
  /// and density difference delta_rho (rho_brine - rho_CO2) [kg/m3].
  Real satNBar(Real h, Real H, Real delta_rho) const;

  const MooseEnum _mode;
  const VariableValue & _sat_n;
  const MaterialProperty<Real> & _H;
  const Real _S_wr;

  // capillary_fringe only
  const PorousFlowCapillaryPressure * _pc_uo;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const Real _gravity_magnitude;

  MaterialProperty<Real> & _h;
};
