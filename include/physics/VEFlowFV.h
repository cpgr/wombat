//* Vertical-equilibrium CCS MOOSE application (wombat)
//*
//* VEFlowFV: cell-centered finite-volume discretization of the VE two-phase flow
//* equations. Creates the MooseVariableFVReal primary variables, the FV flow
//* kernels (VEFVMassTimeDerivative, VEFVAdvectiveFlux), the elemental material
//* chain (VEPorosity, VEPermeability, VESaturation, VEFluidProperties), and the
//* functor materials (VEFVRelPerm, VEFVFluidProperties).
//*
//* Note: no VEGeometry -- the FV kernels build H and grad(z_T) inline from the
//* z_top / z_bottom variables, which the action declares as MooseVariableFVReal.
//* If the user instead pre-declares them, they must be MooseVariableFVReal (a
//* regular FE aux is not reinitialised on FV faces); a mismatched type is caught
//* by MOOSE's variable-type-consistency check. The coupled FV solve also needs
//* SMP full=true, which stays in the user's [Preconditioning] block.
//* Syntax: [Physics/VEFlow/FiniteVolume/<name>].

#pragma once

#include "VEFlowPhysicsBase.h"

class VEFlowFV : public VEFlowPhysicsBase
{
public:
  static InputParameters validParams();

  VEFlowFV(const InputParameters & parameters);

private:
  virtual void addSolverVariables() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addFVKernels() override;
  virtual void addMaterials() override;
  virtual void addFunctorMaterials() override;
  virtual void initializePhysicsAdditional() override;
  virtual InputParameters getAdditionalRMParams() const override;
  virtual void checkGeometryVariableType(const VariableName & var_name) const override;
};
