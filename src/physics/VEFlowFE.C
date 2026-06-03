//* Vertical-equilibrium CCS MOOSE application (wombat)
//* VEFlowFE implementation -- see VEFlowFE.h.

#include "VEFlowFE.h"
#include "MooseVariableBase.h"

registerPhysicsBaseTasks("wombatApp", VEFlowFE);
registerMooseAction("wombatApp", VEFlowFE, "add_variables_physics");
registerMooseAction("wombatApp", VEFlowFE, "add_aux_variable");
registerMooseAction("wombatApp", VEFlowFE, "add_kernel");
registerMooseAction("wombatApp", VEFlowFE, "add_aux_kernel");
registerMooseAction("wombatApp", VEFlowFE, "add_material");

InputParameters
VEFlowFE::validParams()
{
  InputParameters params = VEFlowPhysicsBase::validParams();
  params.addClassDescription(
      "Vertical-equilibrium two-phase (CO2-brine) flow, discretized with the continuous "
      "Galerkin finite element method. Creates the primary variables, the AD flow kernels, "
      "and the FE material chain.");
  params.transferParam<MooseEnum>(MooseVariableBase::validParams(), "order", "variable_order");
  return params;
}

VEFlowFE::VEFlowFE(const InputParameters & parameters)
  : VEFlowPhysicsBase(parameters), _order(getParam<MooseEnum>("variable_order"))
{
}

void
VEFlowFE::addSolverVariables()
{
  const std::string var_type = "MooseVariable";
  for (const auto & var : {_pressure_var, _saturation_var})
  {
    if (!shouldCreateVariable(var, _blocks, /*error_if_aux=*/true))
    {
      reportPotentiallyMissedParameters({"variable_order", "system_names"}, var_type);
      continue;
    }
    auto params = getFactory().getValidParams(var_type);
    assignBlocks(params, _blocks);
    params.set<MooseEnum>("order") = _order;
    params.set<SolverSystemName>("solver_sys") = getSolverSystem(var);
    getProblem().addVariable(var_type, var, params);
  }
}

void
VEFlowFE::addAuxiliaryVariables()
{
  // Geometry elevations as FE LAGRANGE aux variables (a non-zero gradient of z_top is
  // needed for the buoyancy drive). Declared only if the user has not already provided
  // them under these names (e.g. from a mesh field).
  for (const auto & name : {_z_top, _z_bottom})
  {
    if (!shouldCreateVariable(name, _blocks, /*error_if_aux=*/false))
      continue;
    auto params = getFactory().getValidParams("MooseVariable");
    assignBlocks(params, _blocks);
    params.set<MooseEnum>("order") = "FIRST";
    params.set<MooseEnum>("family") = "LAGRANGE";
    getProblem().addAuxVariable("MooseVariable", name, params);
  }

  // Dissolved-CO2 accumulator: elemental (CONSTANT MONOMIAL) because VEDissolvedCO2Aux
  // reads the qp material property ve_dissolution_rate.
  if (_dissolution && shouldCreateVariable(_c_diss, _blocks, /*error_if_aux=*/false))
  {
    auto params = getFactory().getValidParams("MooseVariable");
    assignBlocks(params, _blocks);
    params.set<MooseEnum>("order") = "CONSTANT";
    params.set<MooseEnum>("family") = "MONOMIAL";
    getProblem().addAuxVariable("MooseVariable", _c_diss, params);
  }
}

void
VEFlowFE::addFEKernels()
{
  const auto & g = getParam<RealVectorValue>("gravity");

  // CO2 mass equation (variable = sat_n, phase 0).
  {
    auto params = getFactory().getValidParams("VEMassTimeDerivative");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    params.set<unsigned int>("fluid_phase") = 0;
    getProblem().addKernel("VEMassTimeDerivative", prefix() + "co2_storage", params);
  }
  {
    auto params = getFactory().getValidParams("VEAdvectiveFluxS");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    params.set<unsigned int>("fluid_phase") = 0;
    params.set<std::vector<VariableName>>("pp_top") = {_pressure_var};
    params.set<RealVectorValue>("gravity") = g;
    if (_capillary)
      params.set<bool>("capillary") = true;
    getProblem().addKernel("VEAdvectiveFluxS", prefix() + "co2_flux", params);
  }

  // Brine mass equation (variable = pp_top, phase 1).
  {
    auto params = getFactory().getValidParams("VEMassTimeDerivative");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _pressure_var;
    params.set<unsigned int>("fluid_phase") = 1;
    getProblem().addKernel("VEMassTimeDerivative", prefix() + "brine_storage", params);
  }
  {
    auto params = getFactory().getValidParams("VEAdvectiveFluxP");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _pressure_var;
    params.set<unsigned int>("fluid_phase") = 1;
    params.set<RealVectorValue>("gravity") = g;
    getProblem().addKernel("VEAdvectiveFluxP", prefix() + "brine_flux", params);
  }

  // Optional convective-dissolution sink on the CO2 equation.
  if (_dissolution)
  {
    auto params = getFactory().getValidParams("VEDissolutionSink");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    getProblem().addKernel("VEDissolutionSink", prefix() + "co2_dissolution", params);
  }
}

void
VEFlowFE::addMaterials()
{
  // Geometry: ve_H and ve_grad_z_top from the FE z_top / z_bottom AuxVariables.
  {
    auto params = getFactory().getValidParams("VEGeometry");
    assignBlocks(params, _blocks);
    params.set<std::vector<VariableName>>("z_top") = {_z_top};
    params.set<std::vector<VariableName>>("z_bottom") = {_z_bottom};
    getProblem().addMaterial("VEGeometry", prefix() + "geometry", params);
  }

  addCommonMaterials();

  // Relative permeability (elemental qp property for the FE flux kernels).
  {
    auto params = getFactory().getValidParams("VERelPerm");
    assignBlocks(params, _blocks);
    params.set<UserObjectName>("relperm_uo") = getParam<UserObjectName>("relperm_uo");
    getProblem().addMaterial("VERelPerm", prefix() + "relperm", params);
  }

  // Capillary-pressure chain (capillary = true): VEPlumeReconstruction supplies the
  // sharp-interface plume thickness ve_h that VEUpscaledCapPressure differentiates
  // into the ve_dPcup_dsatn / ve_dPcup_dH coefficients consumed by VEAdvectiveFluxS.
  if (_capillary)
  {
    {
      auto params = getFactory().getValidParams("VEPlumeReconstruction");
      assignBlocks(params, _blocks);
      params.set<std::vector<VariableName>>("sat_n") = {_saturation_var};
      assignSwr(params); // sharp_interface mode is the default
      getProblem().addMaterial("VEPlumeReconstruction", prefix() + "plume_recon", params);
    }
    {
      auto params = getFactory().getValidParams("VEUpscaledCapPressure");
      assignBlocks(params, _blocks);
      params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
      params.set<Real>("pc_entry") = getParam<Real>("pc_entry");
      assignSwr(params);
      getProblem().addMaterial("VEUpscaledCapPressure", prefix() + "cap_pressure", params);
    }
  }
}
