//* Vertical-equilibrium CCS MOOSE application (wombat)
//*
//* VEFlowFE: continuous-Galerkin (finite-element) discretization of the VE
//* two-phase flow equations. Creates the LAGRANGE primary variables, the AD FE
//* flow kernels (VEMassTimeDerivative, VEAdvectiveFluxS/P), and the FE material
//* chain (VEGeometry, VEPorosity, VEPermeability, VESaturation, VEFluidProperties,
//* VERelPerm). Syntax: [Physics/VEFlow/FiniteElement/<name>].

#pragma once

#include "VEFlowPhysicsBase.h"

class VEFlowFE : public VEFlowPhysicsBase
{
public:
  static InputParameters validParams();

  VEFlowFE(const InputParameters & parameters);

private:
  virtual void addSolverVariables() override;
  virtual void addAuxiliaryVariables() override;
  virtual void addFEKernels() override;
  virtual void addMaterials() override;

  /// Variable order/family for the FE primary variables.
  const MooseEnum _order;
};
