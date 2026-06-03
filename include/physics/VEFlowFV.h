//* Vertical-equilibrium CCS MOOSE application (wombat)
//*
//* VEFlowFV: cell-centered finite-volume discretization of the VE two-phase flow
//* equations. Creates the MooseVariableFVReal primary variables, the FV flow
//* kernels (VEFVMassTimeDerivative, VEFVAdvectiveFlux), the elemental material
//* chain (VEPorosity, VEPermeability, VESaturation, VEFluidProperties), and the
//* functor materials (VEFVRelPerm, VEFVFluidProperties).
//*
//* Note: no VEGeometry -- the FV kernels build H and grad(z_T) inline from the
//* z_top / z_bottom FV variables. Those geometry variables MUST themselves be
//* MooseVariableFVReal aux variables (regular FE materials are not reinitialised
//* on FV faces); the Physics only couples them by name. The coupled FV solve also
//* needs SMP full=true, which stays in the user's [Preconditioning] block.
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

  /// Verify z_top / z_bottom are MooseVariableFVReal. Regular FE aux variables are
  /// not reinitialised on FV faces, so they read as zero there and silently kill the
  /// advective flux -- this turns that quiet failure into an up-front error.
  void checkGeometryVariablesAreFV() const;
};
