#pragma once

#include "FunctorMaterial.h"

class SinglePhaseFluidProperties;

/**
 * VEFVFluidProperties
 *
 * FV adapter: publishes the phase density and viscosity as the functors
 * ve_density_n / ve_density_w / ve_viscosity_n / ve_viscosity_w for the FV flux
 * kernel, delegating to two SinglePhaseFluidProperties UserObjects (fp_nw / fp_w)
 * evaluated at the pp_top functor and a constant temperature.
 *
 * The functors evaluate pp_top at whatever face argument the flux kernel supplies.
 * VEFVAdvectiveFlux evaluates them at a central-difference (face-averaged) argument,
 * so the density and viscosity used in the depth-integrated Darcy flux are
 * face-correct -- they depend on the pressure of BOTH adjacent cells and carry AD
 * wrt pp_top on both sides of the face. This is the FV counterpart of the elemental
 * VEFluidProperties material; both wrap the same FP UserObjects at the same pp_top
 * and temperature, so they are consistent by construction.
 *
 * VEFluidProperties is still required in an FV input: the FV mass-storage kernel
 * (VEFVMassTimeDerivative) is elemental and reads the stateful ve_density (with its
 * old value) from VEFluidProperties. This functor adapter feeds only the flux kernel.
 */
class VEFVFluidProperties : public FunctorMaterial
{
public:
  static InputParameters validParams();
  VEFVFluidProperties(const InputParameters & parameters);

protected:
  /// Top-surface pore pressure functor [Pa] -- the EOS evaluation pressure.
  const Moose::Functor<ADReal> & _pp_top;
  /// Constant formation temperature [K], used only to evaluate the EOS.
  const Real _temperature;

  const SinglePhaseFluidProperties & _fp_nw;
  const SinglePhaseFluidProperties & _fp_w;
};
