//* Vertical-equilibrium CCS MOOSE application (wombat)
//* VEFlowFV implementation -- see VEFlowFV.h.

#include "VEFlowFV.h"

registerPhysicsBaseTasks("wombatApp", VEFlowFV);
registerMooseAction("wombatApp", VEFlowFV, "add_variables_physics");
registerMooseAction("wombatApp", VEFlowFV, "add_aux_variable");
registerMooseAction("wombatApp", VEFlowFV, "add_fv_kernel");
registerMooseAction("wombatApp", VEFlowFV, "add_aux_kernel");
registerMooseAction("wombatApp", VEFlowFV, "add_material");
registerMooseAction("wombatApp", VEFlowFV, "add_functor_material");

InputParameters
VEFlowFV::validParams()
{
  InputParameters params = VEFlowPhysicsBase::validParams();
  params.addClassDescription(
      "Vertical-equilibrium two-phase (CO2-brine) flow, discretized with the cell-centered "
      "finite volume method. Creates the FV primary variables, the geometry (and dissolution) "
      "aux variables, the FV flow kernels, the elemental material chain, and the functor "
      "materials. If the user pre-declares z_top / z_bottom they must be MooseVariableFVReal; "
      "SMP full=true stays in the user's [Preconditioning] block.");
  params.addParam<unsigned short>(
      "ghost_layers", 2, "Number of ghosting layers for distributed parallel runs.");

  MooseEnum advected_interp_method("average upwind vanLeer min_mod sou quick venkatakrishnan",
                                   "upwind");
  params.addParam<MooseEnum>(
      "advected_interp_method",
      advected_interp_method,
      "Face interpolation for the advected relative permeability in the FV flux (forwarded to "
      "both VEFVAdvectiveFlux kernels). 'upwind' (default) is first-order; 'average' is "
      "unstabilised central (verification only -- it is also the differentiable choice at "
      "no-flow / zero-velocity states, where the upwind switch is a Jacobian kink); the rest "
      "are second-order TVD limiters that sharpen the plume front (transient only).");

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
VEFlowFV::addAuxiliaryVariables()
{
  // Geometry elevations + dissolved-CO2 accumulator as FV (MooseVariableFVReal) aux
  // variables. Geometry is skipped when define_geometry_variables = false (the user provides
  // z_top / z_bottom, e.g. elemental fields from an Exodus mesh); c_diss is always action-owned
  // (it is an internal accumulator, never a mesh field).
  std::vector<VariableName> aux_names;
  if (_define_geometry_variables)
  {
    checkGeometryNotUserDeclared();
    aux_names.push_back(_z_top);
    aux_names.push_back(_z_bottom);
  }
  if (_dissolution)
    aux_names.push_back(_c_diss);

  for (const auto & name : aux_names)
  {
    auto params = getFactory().getValidParams("MooseVariableFVReal");
    assignBlocks(params, _blocks);
    getProblem().addAuxVariable("MooseVariableFVReal", name, params);
  }
}

void
VEFlowFV::checkGeometryVariableType(const VariableName & var_name) const
{
  // FV needs MooseVariableFVReal geometry: an FE aux is not reinitialised on FV faces, so it
  // reads zero there and silently kills the advective flux.
  if (!isVariableFV(var_name))
    paramError("define_geometry_variables",
               "Geometry variable '",
               var_name,
               "' is not a finite-volume variable. In an FV VEFlow physics z_top/z_bottom must "
               "be MooseVariableFVReal: a regular FE aux variable is not reinitialised on FV "
               "faces, so it reads as zero there and silently kills the advective flux. Declare "
               "it as type = MooseVariableFVReal, or set define_geometry_variables=true to let "
               "the action declare it.");
}

void
VEFlowFV::addFVKernels()
{
  const auto & g = getParam<RealVectorValue>("gravity");
  const VariableName & z_top = _z_top;
  const VariableName & z_bottom = _z_bottom;

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
    params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("advected_interp_method");
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
    params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("advected_interp_method");
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
      params.set<MooseFunctorName>("z_top") = _z_top;
      params.set<MooseFunctorName>("z_bottom") = _z_bottom;
      params.set<Real>("S_wr") = getParam<Real>("S_wr");
      params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
    }
    getProblem().addFunctorMaterial("VEFVFluidProperties", prefix() + "fv_fluid_props", params);
  }

  // Capillary-pressure functors (capillary != none). delta_rho = rho_w - rho_n is read
  // from the ve_density_n / ve_density_w functors above (the defaults of density_nw /
  // density_w), so no separate density input is needed; VEFVCapPressure computes h
  // internally (sharp closed form, or capillary_fringe Newton inversion with pc_uo), so
  // -- unlike the FE side -- no VEPlumeReconstruction is required.
  if (_capillary)
  {
    auto params = getFactory().getValidParams("VEFVCapPressure");
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>("sat_n") = _saturation_var;
    params.set<MooseFunctorName>("z_top") = _z_top;
    params.set<MooseFunctorName>("z_bottom") = _z_bottom;
    params.set<Real>("pc_entry") = getParam<Real>("pc_entry");
    params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
    assignSwr(params);
    assignCapillaryMode(params);
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
