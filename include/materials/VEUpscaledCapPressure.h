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
 * CO2 column of height h (sharp-interface upscaled Pc; Nordbotten & Celia). The
 * optional constant pc_entry adds a capillary entry/fringe offset at the plume
 * top -- a minimal capillary-fringe correction. A full Nordbotten & Dahle fringe
 * integral (consuming the column Pc UserObject, paralleling the capillary_fringe
 * mode of VEPlumeReconstruction) is a future extension.
 *
 * Reads ve_h (from VEPlumeReconstruction) and ve_density. Like ve_h, ve_pc_up is
 * currently a plain (non-AD) Real diagnostic property: it is not yet in the
 * nonlinear residual. When ve_h is upgraded to an ADMaterialProperty (see the
 * VEPlumeReconstruction notes) so that Pc^up can enter a two-pressure flux, this
 * property should likewise become AD so d(Pc^up)/d(sat_n) is captured.
 */
class VEUpscaledCapPressure : public Material
{
public:
  static InputParameters validParams();
  VEUpscaledCapPressure(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Scalar magnitude of gravity [m/s2].
  const Real _gravity_magnitude;
  /// Capillary entry/fringe pressure added at the plume top [Pa].
  const Real _pc_entry;

  /// Reconstructed CO2 plume thickness h [m] (from VEPlumeReconstruction).
  const MaterialProperty<Real> & _h;
  /// Depth-averaged phase densities [kg/m3]: [0] = CO2, [1] = brine.
  const ADMaterialProperty<std::vector<Real>> & _density;

  /// Upscaled capillary pressure Pc^up [Pa].
  MaterialProperty<Real> & _pc_up;
};
