//* Vertical-equilibrium CCS MOOSE application (wombat)
//* VEFlowFV implementation -- see VEFlowFV.h.

#include "VEFlowFV.h"

registerPhysicsBaseTasks("wombatApp", VEFlowFV);
registerMooseAction("wombatApp", VEFlowFV, "add_variables_physics");
registerMooseAction("wombatApp", VEFlowFV, "add_fv_kernel");
registerMooseAction("wombatApp", VEFlowFV, "add_material");
registerMooseAction("wombatApp", VEFlowFV, "add_functor_material");

InputParameters
VEFlowFV::validParams()
{
  InputParameters params = VEFlowPhysicsBase::validParams();
  params.addClassDescription(
      "Vertical-equilibrium two-phase (CO2-brine) flow, discretized with the cell-centered "
      "finite volume method. Creates the FV primary variables, the FV flow kernels, the "
      "elemental material chain, and the functor materials. z_top / z_bottom must be FV "
      "aux variables; SMP full=true stays in the user's [Preconditioning] block.");
  params.addParam<unsigned short>(
      "ghost_layers", 2, "Number of ghosting layers for distributed parallel runs.");
  params.addParamNamesToGroup("ghost_layers", "Advanced");
  return params;
}

VEFlowFV::VEFlowFV(const InputParameters & parameters) : VEFlowPhysicsBase(parameters) {}

void
VEFlowFV::initializePhysicsAdditional()
{
  getProblem().needFV();
}

void
VEFlowFV::addSolverVariables()
{
  const std::string var_type = "MooseVariableFVReal";
  for (const auto & var : {_pressure_var, _saturation_var})
  {
    if (!shouldCreateVariable(var, _blocks, /*error_if_aux=*/true))
    {
      reportPotentiallyMissedParameters({"system_names"}, var_type);
      continue;
    }
    auto params = getFactory().getValidParams(var_type);
    assignBlocks(params, _blocks);
    params.set<SolverSystemName>("solver_sys") = getSolverSystem(var);
    getProblem().addVariable(var_type, var, params);
  }
}

void
VEFlowFV::checkGeometryVariablesAreFV() const
{
  for (const std::string & p : {std::string("z_top"), std::string("z_bottom")})
  {
    const auto names = getParam<std::vector<VariableName>>(p);
    if (names.empty())
      paramError(p,
                 "In the finite-volume VE flow physics, '",
                 p,
                 "' must couple to a MooseVariableFVReal AuxVariable (a constant is not "
                 "allowed): the FV kernels difference its face values to build H and grad(z_T).");

    const auto & name = names[0];
    if (!getProblem().hasVariable(name))
      paramError(p, "Coupled geometry variable '", name, "' does not exist.");
    if (!isVariableFV(name))
      paramError(p,
                 "Coupled geometry variable '",
                 name,
                 "' is not a finite-volume variable. In the FV VE flow physics, z_top and "
                 "z_bottom must be MooseVariableFVReal: regular FE aux variables are not "
                 "reinitialised on FV faces, so they read as zero at every face and silently "
                 "kill the advective flux. Declare them with type = MooseVariableFVReal.");
  }
}

void
VEFlowFV::addFVKernels()
{
  checkGeometryVariablesAreFV();

  const auto & g = getParam<RealVectorValue>("gravity");
  const VariableName z_top = getParam<std::vector<VariableName>>("z_top")[0];
  const VariableName z_bottom = getParam<std::vector<VariableName>>("z_bottom")[0];

  // CO2 mass equation (variable = sat_n, phase 0).
  {
    auto params = getFactory().getValidParams("VEFVMassTimeDerivative");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    params.set<unsigned int>("fluid_phase") = 0;
    params.set<std::vector<VariableName>>("z_top") = {z_top};
    params.set<std::vector<VariableName>>("z_bottom") = {z_bottom};
    getProblem().addFVKernel("VEFVMassTimeDerivative", prefix() + "co2_storage", params);
  }
  {
    auto params = getFactory().getValidParams("VEFVAdvectiveFlux");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    params.set<unsigned int>("fluid_phase") = 0;
    params.set<MooseFunctorName>("pp_top") = _pressure_var;
    params.set<MooseFunctorName>("z_top") = z_top;
    params.set<MooseFunctorName>("z_bottom") = z_bottom;
    params.set<RealVectorValue>("gravity") = g;
    if (_capillary)
      params.set<bool>("capillary") = true;
    getProblem().addFVKernel("VEFVAdvectiveFlux", prefix() + "co2_flux", params);
  }

  // Brine mass equation (variable = pp_top, phase 1).
  {
    auto params = getFactory().getValidParams("VEFVMassTimeDerivative");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _pressure_var;
    params.set<unsigned int>("fluid_phase") = 1;
    params.set<std::vector<VariableName>>("z_top") = {z_top};
    params.set<std::vector<VariableName>>("z_bottom") = {z_bottom};
    getProblem().addFVKernel("VEFVMassTimeDerivative", prefix() + "brine_storage", params);
  }
  {
    auto params = getFactory().getValidParams("VEFVAdvectiveFlux");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _pressure_var;
    params.set<unsigned int>("fluid_phase") = 1;
    params.set<MooseFunctorName>("pp_top") = _pressure_var;
    params.set<MooseFunctorName>("z_top") = z_top;
    params.set<MooseFunctorName>("z_bottom") = z_bottom;
    params.set<RealVectorValue>("gravity") = g;
    getProblem().addFVKernel("VEFVAdvectiveFlux", prefix() + "brine_flux", params);
  }

  // Optional convective-dissolution sink on the CO2 equation.
  if (_dissolution)
  {
    auto params = getFactory().getValidParams("VEFVDissolutionSink");
    assignBlocks(params, _blocks);
    params.set<NonlinearVariableName>("variable") = _saturation_var;
    getProblem().addFVKernel("VEFVDissolutionSink", prefix() + "co2_dissolution", params);
  }
}

void
VEFlowFV::addMaterials()
{
  // No VEGeometry in FV: the kernels build H and grad(z_T) inline from z_top/z_bottom.
  addCommonMaterials();
}

void
VEFlowFV::addFunctorMaterials()
{
  // Relative permeability functors, evaluated at the upwind face.
  {
    auto params = getFactory().getValidParams("VEFVRelPerm");
    assignBlocks(params, _blocks);
    params.set<UserObjectName>("relperm_uo") = getParam<UserObjectName>("relperm_uo");
    params.set<MooseFunctorName>("sat_n") = _saturation_var;
    getProblem().addFunctorMaterial("VEFVRelPerm", prefix() + "relperm", params);
  }

  // Face-correct fluid-property functors for the FV advective flux.
  {
    auto params = getFactory().getValidParams("VEFVFluidProperties");
    assignBlocks(params, _blocks);
    params.set<UserObjectName>("fp_nw") = getParam<UserObjectName>("fp_nw");
    params.set<UserObjectName>("fp_w") = getParam<UserObjectName>("fp_w");
    params.set<MooseFunctorName>("pp_top") = _pressure_var;
    params.set<Real>("temperature") = getParam<Real>("temperature");
    params.set<MooseEnum>("eos_reference_depth") = getParam<MooseEnum>("eos_reference_depth");
    if (_interface_eos)
    {
      params.set<MooseFunctorName>("sat_n") = _saturation_var;
      params.set<MooseFunctorName>("z_top") = getParam<std::vector<VariableName>>("z_top")[0];
      params.set<MooseFunctorName>("z_bottom") = getParam<std::vector<VariableName>>("z_bottom")[0];
      params.set<Real>("S_wr") = getParam<Real>("S_wr");
      params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
    }
    getProblem().addFunctorMaterial("VEFVFluidProperties", prefix() + "fv_fluid_props", params);
  }

  // Capillary-pressure functors (capillary = true). delta_rho = rho_w - rho_n is read
  // from the ve_density_n / ve_density_w functors above (the defaults of density_nw /
  // density_w), so no separate density input is needed; VEFVCapPressure computes the
  // sharp-interface h internally, so -- unlike the FE side -- no VEPlumeReconstruction
  // is required.
  if (_capillary)
  {
    auto params = getFactory().getValidParams("VEFVCapPressure");
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>("sat_n") = _saturation_var;
    params.set<MooseFunctorName>("z_top") = getParam<std::vector<VariableName>>("z_top")[0];
    params.set<MooseFunctorName>("z_bottom") = getParam<std::vector<VariableName>>("z_bottom")[0];
    params.set<Real>("pc_entry") = getParam<Real>("pc_entry");
    params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
    assignSwr(params);
    getProblem().addFunctorMaterial("VEFVCapPressure", prefix() + "cap_pressure", params);
  }
}

InputParameters
VEFlowFV::getAdditionalRMParams() const
{
  // Match the advective-flux kernel's ghosting so functor gradient reconstruction
  // and upwind face evaluation work in parallel.
  const auto necessary_layers =
      std::max(getParam<unsigned short>("ghost_layers"), (unsigned short)2);
  InputParameters params = getFactory().getValidParams("VEFVAdvectiveFlux");
  params.set<unsigned short>("ghost_layers") = necessary_layers;
  return params;
}
