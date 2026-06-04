//* Vertical-equilibrium CCS MOOSE application (wombat)
//*
//* VEFlowPhysicsBase: common parameters and shared object-creation helpers for
//* the VE two-phase (CO2-brine) flow Physics. Concrete discretizations are
//* VEFlowFE (continuous Galerkin) and VEFlowFV (cell-centered finite volume),
//* mirroring the framework's DiffusionPhysicsBase / DiffusionCG / DiffusionFV
//* pattern. The Physics adds the two primary variables (pp_top, sat_n), the four
//* depth-integrated flow kernels (CO2 + brine, storage + flux), and the standard
//* material chain, so a standard VE input no longer hand-writes [Variables],
//* [Kernels]/[FVKernels], [Materials], or [FunctorMaterials].

#pragma once

#include "PhysicsBase.h"

/**
 * Base class hosting the parameters and shared helpers common to the VE flow
 * Physics. Abstract: the variable type and the kernel/material wiring that
 * differ between FE and FV are supplied by the derived classes.
 */
class VEFlowPhysicsBase : public PhysicsBase
{
public:
  static InputParameters validParams();

  VEFlowPhysicsBase(const InputParameters & parameters);

protected:
  /// Name of the pore-pressure primary variable (brine mass equation).
  const VariableName & _pressure_var;
  /// Name of the depth-averaged CO2 saturation primary variable (CO2 mass equation).
  const VariableName & _saturation_var;
  /// Name of the top-surface elevation z_T aux variable (declared by the action).
  const VariableName & _z_top;
  /// Name of the bottom-surface elevation z_B aux variable (declared by the action).
  const VariableName & _z_bottom;
  /// Name of the dissolved-CO2 areal-mass aux variable (declared when dissolution is on).
  const VariableName & _c_diss;
  /// True when the action declares z_top/z_bottom; false when the user provides them.
  const bool _define_geometry_variables;
  /// True when eos_reference_depth = interface (rho/mu at the CO2-brine contact).
  const bool _interface_eos;
  /// True when the convective-dissolution chain (material + aux + sink) should be added.
  const bool _dissolution;
  /// True when the capillary-pressure term and its material(s) should be added.
  const bool _capillary;

  /// Forward a coupled-var-or-constant Physics parameter to a downstream object.
  void assignCoupled(InputParameters & params, const std::string & name) const;

  /// Forward S_wr to a downstream object only if the user set it (else the object
  /// keeps its own default). Shared by the interface-EOS and capillary materials.
  void assignSwr(InputParameters & params) const;

  /// Add the elemental materials shared by FE and FV: VEPorosity, VEPermeability,
  /// VESaturation, the elemental VEFluidProperties, and (when dissolution is on) the
  /// VEDissolution rate material.
  void addCommonMaterials();

  /// Add the VEDissolvedCO2Aux accumulator when dissolution is on (shared FE/FV).
  virtual void addAuxiliaryKernels() override;

  /// Error (with an actionable "rebuild the Exodus file" message) if any field variable the
  /// physics consumes -- name-valued phi_bar / K_up_*, the geometry variables, c_diss -- is
  /// absent. Called at add_material, once all variables have been created.
  void checkRequiredFields();

  /// Verify a present geometry variable (z_top/z_bottom) is the right kind for this
  /// discretization, with an actionable paramError if not. FE requires a nodal variable
  /// (for the buoyancy gradient); FV requires a MooseVariableFVReal.
  virtual void checkGeometryVariableType(const VariableName & var_name) const = 0;

private:
  /// Add the elemental VEFluidProperties material (ve_density / ve_viscosity).
  void addFluidPropertiesMaterial();
};
