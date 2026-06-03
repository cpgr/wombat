//* Vertical-equilibrium CCS MOOSE application (wombat)
//* VEFlowPhysicsBase implementation -- see VEFlowPhysicsBase.h.

#include "VEFlowPhysicsBase.h"

InputParameters
VEFlowPhysicsBase::validParams()
{
  InputParameters params = PhysicsBase::validParams();
  params.addClassDescription(
      "Base class for the vertical-equilibrium two-phase (CO2-brine) flow Physics. "
      "Adds the primary variables, the four depth-integrated flow kernels, and the "
      "standard VE material chain.");

  // Primary variable names (overridable so several Physics can coexist).
  params.addParam<VariableName>(
      "pressure_variable", "pp_top", "Name of the pore-pressure-at-top primary variable.");
  params.addParam<VariableName>(
      "saturation_variable", "sat_n", "Name of the depth-averaged CO2 saturation primary variable.");

  // Geometry (coupled AuxVariables, typically nodal/cell fields from the upscaling mesh).
  params.addRequiredCoupledVar("z_top", "Top-surface elevation z_T [m].");
  params.addRequiredCoupledVar("z_bottom", "Bottom-surface elevation z_B [m].");

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

  params.addParam<bool>(
      "enable_dissolution",
      false,
      "If true, add the convective-dissolution sink on the CO2 mass equation. Requires a "
      "VEDissolution material in [Materials] providing ve_dissolution_rate.");

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
    _interface_eos(getParam<MooseEnum>("eos_reference_depth") == "interface"),
    _dissolution(getParam<bool>("enable_dissolution")),
    _capillary(getParam<bool>("capillary"))
{
  saveSolverVariableName(_pressure_var);
  saveSolverVariableName(_saturation_var);

  if (_interface_eos && !isParamValid("S_wr"))
    paramError("S_wr", "S_wr is required when eos_reference_depth = interface.");
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
VEFlowPhysicsBase::addCommonMaterials()
{
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
    assignCoupled(params, "z_top");
    assignCoupled(params, "z_bottom");
    params.set<Real>("S_wr") = getParam<Real>("S_wr");
    params.set<RealVectorValue>("gravity") = getParam<RealVectorValue>("gravity");
  }

  getProblem().addMaterial("VEFluidProperties", prefix() + "fluid_props", params);
}
