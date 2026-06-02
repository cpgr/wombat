# Jacobian verification of the FV advective flux with a pressure-dependent EOS.
#
# This is the variable-density companion to the constant-property FV tests: it
# exercises the face-correct density path that those tests cannot, because for
# constant fluid properties the face average of density reduces to a constant and
# the rework is a no-op. Here the non-wetting phase is real CO2 (CO2FluidProperties)
# and the wetting phase a fixed-salinity brine (SimpleBrineFluidProperties), so both
# densities and viscosities vary with pressure across every face. A pressure gradient
# and a sloped top surface make the mass-flux coefficient AND the buoyancy term
# density-dependent.
#
# VEFVFluidProperties publishes density/viscosity as functors; VEFVAdvectiveFlux
# evaluates them at the central-difference (face-averaged) argument, which carries AD
# wrt pp_top on both sides of the face. PetscJacobianTester confirms the hand-assembled
# AD Jacobian matches the finite-difference reference -- i.e. d(rho)/d(pp_top) and
# d(mu)/d(pp_top) are seeded correctly through the face evaluation.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  VEDictator = ve_dictator
[]

[FluidProperties]
  [co2_fp]
    type = CO2FluidProperties
  []
  [brine_fp]
    type = SimpleBrineFluidProperties
    salt_mass_fraction = 0.1
  []
[]

[UserObjects]
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = sharp_interface
  []
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.1
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '100.0 - 0.1 * x' # sloped top -> nonzero grad(z_T) buoyancy drive
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '80.0 - 0.1 * x'
  []
  [pp_top_func]
    type = ParsedFunction
    expression = '12.0e6 - 4.0e4 * x' # 12 MPa -> 8 MPa: supercritical CO2, rho varies
  []
[]

[Variables]
  [pp_top]
    type = MooseVariableFVReal
  []
  [sat_n]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_func
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5 # off the sharp-relperm kink, so CO2 mobility and its derivative are live
  []
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = z_top_func
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_func
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
  []
[]

[FVKernels]
  [co2_storage]
    type = VEFVMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
    z_top = z_top
    z_bottom = z_bottom
  []
  [co2_flux]
    type = VEFVAdvectiveFlux
    variable = sat_n
    fluid_phase = 0
    pp_top = pp_top
    z_top = z_top
    z_bottom = z_bottom
  []
  [brine_storage]
    type = VEFVMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
    z_top = z_top
    z_bottom = z_bottom
  []
  [brine_flux]
    type = VEFVAdvectiveFlux
    variable = pp_top
    fluid_phase = 1
    pp_top = pp_top
    z_top = z_top
    z_bottom = z_bottom
  []
[]

[FVBCs]
  [pp_left]
    type = FVFunctionDirichletBC
    variable = pp_top
    boundary = left
    function = pp_top_func
  []
  [pp_right]
    type = FVFunctionDirichletBC
    variable = pp_top
    boundary = right
    function = pp_top_func
  []
[]

[Materials]
  [porosity]
    type = VEPorosity
    phi_bar = 0.2
  []
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
  []
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
    temperature = 323.15
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
[]

[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n = sat_n
  []
  [fv_fluid_props]
    type = VEFVFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
    temperature = 323.15
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 1.0
  dt = 1.0
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]
