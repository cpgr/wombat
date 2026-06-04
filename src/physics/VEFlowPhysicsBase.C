//* Vertical-equilibrium CCS MOOSE application (wombat)
//* VEFlowPhysicsBase implementation -- see VEFlowPhysicsBase.h.

#include "VEFlowPhysicsBase.h"
#include "MooseUtils.h"
#include "AddAuxVariableAction.h"

InputParameters
VEFlowPhysicsBase::validParams()
{
  InputParameters params = PhysicsBase::validParams();
  params.addClassDescription(
      "Base class for the vertical-equilibrium two-phase (CO2-brine) flow Physics. "
      "Adds the two primary variables, the four depth-integrated flow kernels, the standard "
      "VE material chain, and the geometry / dissolution auxiliary variables. The user "
      "references the PRIMARY VARIABLE NAMES -- 'pp_top' (pressure) and 'sat_n' (CO2 "
      "saturation) by default -- when writing BCs, injection wells, ICs and postprocessors; "
      "rename them with pressure_variable / saturation_variable if needed.");

  // Primary variable names. NOTE: these are the names the user must use for BCs, wells,
  // ICs and postprocessors -- the action only creates the variables, it does not own the
  // boundary/well/IC setup (which is model-specific).
  params.addParam<VariableName>(
      "pressure_variable", "pp_top", "Name of the pore-pressure-at-top primary variable.");
  params.addParam<VariableName>(
      "saturation_variable", "sat_n", "Name of the depth-averaged CO2 saturation primary variable.");

  // Geometry. The action DECLARES these aux variables (with the correct type for the
  // discretization) under the names below if they do not already exist; the user supplies
  // their values via an IC or by reading them from the mesh. Override a name to point at an
  // existing field (e.g. an Exodus nodal field), in which case the action skips declaration.
  params.addParam<VariableName>(
      "z_top", "z_top", "Name of the top-surface elevation z_T [m] aux variable.");
  params.addParam<VariableName>(
      "z_bottom", "z_bottom", "Name of the bottom-surface elevation z_B [m] aux variable.");
  params.addParam<bool>(
      "define_geometry_variables",
      true,
      "If true (default) the action declares the z_top / z_bottom aux variables with the "
      "correct type for the discretization (FE LAGRANGE / FV MooseVariableFVReal); you supply "
      "their values via an IC. If false the action does NOT declare them -- you provide them "
      "yourself (e.g. read elemental fields from an Exodus mesh with initial_from_file_var) -- "
      "and the action checks that they exist and are the right type, erroring otherwise. Set "
      "false whenever the geometry comes from the mesh, so you control the variables explicitly "
      "rather than relying on the action's declaration merging with yours.");

  // Petrophysics (constant for verification, or coupled fields on real cases).
  params.addRequiredCoupledVar("phi_bar", "Depth-averaged porosity [-].");
  params.addRequiredCoupledVar("K_up_xx", "Depth-integrated permeability K_up xx-component [m3].");
  params.addRequiredCoupledVar("K_up_yy", "Depth-integrated permeability K_up yy-component [m3].");
  params.addCoupledVar("K_up_xy", "Depth-integrated permeability K_up xy-component [m3] (default 0).");

  // Fluids.
  params.addRequiredParam<UserObjectName>(
      "fp_nw", "SinglePhaseFluidProperties UserObject for the non-wetting phase (CO2).");
  params.addRequiredParam<UserObjectName>(
      "fp_w", "SinglePhaseFluidProperties UserObject for the wetting phase (brine).");
  params.addRequiredParam<UserObjectName>(
      "relperm_uo", "VERelativePermeability UserObject providing the kr_c(sat_n) curves.");
  params.addRangeCheckedParam<Real>(
      "temperature", 313.15, "temperature > 0", "Isothermal formation temperature [K].");

  MooseEnum eos_reference_depth("top_surface interface", "top_surface");
  params.addParam<MooseEnum>(
      "eos_reference_depth",
      eos_reference_depth,
      "Depth at which the EOS pressure is evaluated (CLAUDE.md key-subtlety 2). top_surface "
      "(default): rho/mu at pp_top. interface: at the CO2-brine contact z_T + h; requires S_wr.");
  params.addRangeCheckedParam<Real>(
      "S_wr",
      "S_wr >= 0 & S_wr < 1",
      "Residual water saturation in the CO2 zone [-]. Required for eos_reference_depth=interface.");

  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(0.0, 0.0, -9.81),
      "Gravity vector [m/s2]. Only the magnitude enters the buoyancy and interface-EOS terms.");

  // Convective dissolution. enable_dissolution = true fully wires it: the action creates
  // the VEDissolution material, the c_diss dissolved-mass aux variable, the VEDissolvedCO2Aux
  // accumulator, and the VE*DissolutionSink kernel. dissolution_flux is then REQUIRED.
  params.addParam<bool>(
      "enable_dissolution",
      false,
      "If true, create the full convective-dissolution chain: the VEDissolution material, the "
      "dissolved-CO2 aux variable + VEDissolvedCO2Aux accumulator, and the sink kernel. "
      "dissolution_flux is required when this is true.");
  params.addParam<VariableName>(
      "dissolved_co2_variable",
      "c_diss",
      "Name of the areal dissolved-CO2 mass aux variable [kg/m2] declared when "
      "enable_dissolution = true.");
  params.addRangeCheckedParam<Real>(
      "dissolution_flux",
      "dissolution_flux >= 0",
      "Constant-flux convective dissolution rate q0 [kg/m2/s] (VEDissolution). REQUIRED when "
      "enable_dissolution = true.");
  params.addRangeCheckedParam<Real>(
      "dissolution_s_ref",
      "dissolution_s_ref > 0",
      "Gate reference CO2 saturation [-] (VEDissolution s_ref). Optional; defaults to the "
      "material default (0.05) if unset.");
  params.addRangeCheckedParam<Real>(
      "dissolution_c_cap",
      "dissolution_c_cap > 0",
      "Optional column CO2 capacity [kg/m2] (VEDissolution c_cap); when set, the action wires "
      "the lagged dissolved_co2 = c_diss so dissolution tapers to zero as c_diss -> c_cap.");
  params.addParamNamesToGroup("enable_dissolution dissolved_co2_variable dissolution_flux "
                              "dissolution_s_ref dissolution_c_cap",
                              "Dissolution");

  params.addParam<bool>(
      "capillary",
      false,
      "If true, add the upscaled capillary-pressure term grad(Pc^up) to the CO2 Darcy "
      "potential (two-pressure VE) and create the cap-pressure material(s) that supply it. "
      "Sharp interface: Pc^up = (rho_w - rho_n)*|g|*sat_n*H/(1 - S_wr) + pc_entry.");
  params.addParam<Real>(
      "pc_entry", 0.0, "Capillary entry/fringe pressure offset [Pa] (used when capillary = true).");

  return params;
}

VEFlowPhysicsBase::VEFlowPhysicsBase(const InputParameters & parameters)
  : PhysicsBase(parameters),
    _pressure_var(getParam<VariableName>("pressure_variable")),
    _saturation_var(getParam<VariableName>("saturation_variable")),
    _z_top(getParam<VariableName>("z_top")),
    _z_bottom(getParam<VariableName>("z_bottom")),
    _c_diss(getParam<VariableName>("dissolved_co2_variable")),
    _define_geometry_variables(getParam<bool>("define_geometry_variables")),
    _interface_eos(getParam<MooseEnum>("eos_reference_depth") == "interface"),
    _dissolution(getParam<bool>("enable_dissolution")),
    _capillary(getParam<bool>("capillary"))
{
  saveSolverVariableName(_pressure_var);
  saveSolverVariableName(_saturation_var);

  if (_interface_eos && !isParamValid("S_wr"))
    paramError("S_wr", "S_wr is required when eos_reference_depth = interface.");

  if (_dissolution && !isParamValid("dissolution_flux"))
    paramError("dissolution_flux",
               "dissolution_flux is required when enable_dissolution = true.");
}

void
VEFlowPhysicsBase::assignCoupled(InputParameters & params, const std::string & name) const
{
  // Forward a coupled-var-or-constant value to the downstream object's coupled-var
  // parameter, using MOOSE's own parameter-transfer machinery. This copies both the
  // AuxVariable-name value (and marks it valid) and any constant "default coupled
  // value" (e.g. phi_bar = 0.2). Skip when the Physics has no value for this name so
  // an unset optional (e.g. K_up_xy) keeps the downstream object's own default.
  const auto & src = parameters();
  if (!src.isParamValid(name) && !src.hasDefaultCoupledValue(name))
    return;
  params.applySpecificParameters(src, {name}, /*allow_private=*/true);
}

void
VEFlowPhysicsBase::assignSwr(InputParameters & params) const
{
  if (isParamValid("S_wr"))
    params.set<Real>("S_wr") = getParam<Real>("S_wr");
}

void
VEFlowPhysicsBase::checkRequiredFields()
{
  // Validate that every field variable this physics consumes actually exists. Catches the
  // common real-field mistake of pointing phi_bar / K_up / geometry at an Exodus field name
  // that the upscaling workflow did not emit (or that is misspelled), with a message that
  // tells the user to rebuild the mesh -- instead of a downstream "no such variable" or a
  // silently-zero field. Runs at add_material, after all variables (action-declared,
  // user [AuxVariables], and initial_from_file_var) have been created. Coupled params given
  // as a numeric constant (e.g. phi_bar = 0.2) carry no variable name and are skipped.
  std::vector<std::string> missing;
  auto check_coupled = [&](const std::string & pname)
  {
    // Skip an unset optional coupled var (e.g. K_up_xy) or one given as a constant
    // (no variable name to resolve); getParam would throw on a never-set coupled var.
    if (!isParamValid(pname))
      return;
    for (const auto & vn : getParam<std::vector<VariableName>>(pname))
      if (!MooseUtils::parsesToReal(vn) && !getProblem().hasVariable(vn))
        missing.push_back(std::string(vn) + " (" + pname + ")");
  };
  check_coupled("phi_bar");
  check_coupled("K_up_xx");
  check_coupled("K_up_yy");
  check_coupled("K_up_xy");

  for (const auto & vn : {_z_top, _z_bottom})
    if (!getProblem().hasVariable(vn))
      missing.push_back(std::string(vn) + " (geometry)");
  if (_dissolution && !getProblem().hasVariable(_c_diss))
    missing.push_back(std::string(_c_diss) + " (dissolved_co2_variable)");

  // Type-correctness of the geometry variables for this discretization (FV must be
  // MooseVariableFVReal, FE must be nodal). Meaningful chiefly when
  // define_geometry_variables = false (the user supplies them); when the action declares
  // them it declares the right type, so this passes. Only checked when present -- a missing
  // one is reported below.
  for (const auto & vn : {_z_top, _z_bottom})
    if (getProblem().hasVariable(vn))
      checkGeometryVariableType(vn);

  if (missing.empty())
    return;

  std::string list;
  for (const auto & m : missing)
    list += (list.empty() ? "" : ", ") + m;
  mooseError(
      "VEFlow: the following field variable(s) consumed by this physics are not present in "
      "the problem: ",
      list,
      ".\nIf your model reads petrophysical/geometry fields from an Exodus mesh, rebuild the "
      "mesh file so it contains element or nodal fields with these exact names -- the upscaling "
      "workflow must emit z_top, z_bottom, phi_bar, K_up_xx, K_up_yy (and K_up_xy if anisotropic) "
      "-- and read them via [AuxVariables] with initial_from_file_var. Otherwise supply them as "
      "constants or declare the variables yourself.");
}

void
VEFlowPhysicsBase::addCommonMaterials()
{
  // Fail early with an actionable message if any consumed field is missing.
  checkRequiredFields();

  // Porosity.
  {
    auto params = getFactory().getValidParams("VEPorosity");
    assignBlocks(params, _blocks);
    assignCoupled(params, "phi_bar");
    getProblem().addMaterial("VEPorosity", prefix() + "porosity", params);
  }

  // Permeability tensor.
  {
    auto params = getFactory().getValidParams("VEPermeability");
    assignBlocks(params, _blocks);
    assignCoupled(params, "K_up_xx");
    assignCoupled(params, "K_up_yy");
    assignCoupled(params, "K_up_xy"); // no-op if the user left it at the default
    getProblem().addMaterial("VEPermeability", prefix() + "permeability", params);
  }

  // Saturation vector [sat_n, 1 - sat_n].
  {
    auto params = getFactory().getValidParams("VESaturation");
    assignBlocks(params, _blocks);
    params.set<std::vector<VariableName>>("sat_n") = {_saturation_var};
    getProblem().addMaterial("VESaturation", prefix() + "saturation", params);
  }

  // Elemental fluid properties (needed by the mass-storage kernels in FE and FV).
  addFluidPropertiesMaterial();

  // Convective-dissolution rate material (shared by FE and FV).
  if (_dissolution)
  {
    auto params = getFactory().getValidParams("VEDissolution");
    assignBlocks(params, _blocks);
    params.set<Real>("dissolution_flux") = getParam<Real>("dissolution_flux");
    if (isParamValid("dissolution_s_ref"))
      params.set<Real>("s_ref") = getParam<Real>("dissolution_s_ref");
    if (isParamValid("dissolution_c_cap"))
    {
      params.set<Real>("c_cap") = getParam<Real>("dissolution_c_cap");
      params.set<std::vector<VariableName>>("dissolved_co2") = {_c_diss};
    }
    getProblem().addMaterial("VEDissolution", prefix() + "dissolution", params);
  }
}

void
VEFlowPhysicsBase::checkGeometryNotUserDeclared() const
{
  // With define_geometry_variables = true the action owns z_top/z_bottom; a user [AuxVariables]
  // declaration of the same name is a conflict (it would otherwise silently merge if the types
  // match, or trip an opaque "Mismatching types" error if not). Detect it reliably from the
  // ActionWarehouse (populated at parse time, independent of action execution order) and give a
  // clear message pointing at the fix.
  const auto & aux_actions = _awh.getActions<AddAuxVariableAction>();
  for (const auto & vn : {_z_top, _z_bottom})
    for (const auto * aux_action : aux_actions)
      if (aux_action->name() == vn)
        paramError("define_geometry_variables",
                   "AuxVariable '",
                   vn,
                   "' has already been declared in [AuxVariables], but this physics also "
                   "declares it because define_geometry_variables = true. Remove the "
                   "[AuxVariables/",
                   vn,
                   "] block, or set define_geometry_variables = false to provide the geometry "
                   "variables yourself.");
}

void
VEFlowPhysicsBase::addAuxiliaryKernels()
{
  // Accumulate the areal dissolved-CO2 mass when dissolution is active. The c_diss aux
  // variable is declared by the derived class (correct type per discretization).
  if (!_dissolution)
    return;
  auto params = getFactory().getValidParams("VEDissolvedCO2Aux");
  assignBlocks(params, _blocks);
  params.set<AuxVariableName>("variable") = _c_diss;
  getProblem().addAuxKernel("VEDissolvedCO2Aux", prefix() + "dissolved_co2", params);
}

void
VEFlowPhysicsBase::addFluidPropertiesMaterial()
{
  auto params = getFactory().getValidParams("VEFluidProperties");
  assignBlocks(params, _blocks);
  params.set<UserObjectName>("fp_nw") = getParam<UserObjectName>("fp_nw");
  params.set<UserObjectName>("fp_w") = getParam<UserObjectName>("fp_w");
  params.set<std::vector<VariableName>>("pp_top") = {_pressure_var};
  params.set<Real>("temperature") = getParam<Real>("temperature");
  params.set<MooseEnum>("eos_reference_depth") = getParam<MooseEnum>("eos_reference_depth");

  if (_interface_eos)
  {
    params.set<std::vector<VariableName>>("sat_n") = {_saturation_var};
    params.set<std::vector<VariableName>>("z_top") = {_z_top};
    params.set<std::vector<VariableName>>("z_bottom") = {_z_bottom};
    params.set<Real>("S_wr") = getParam<Real>("S_wr");
    params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
  }

  getProblem().addMaterial("VEFluidProperties", prefix() + "fluid_props", params);
}
