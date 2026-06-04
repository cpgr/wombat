# VE flow [Physics] -- FV, define_geometry_variables = false (user owns z_top/z_bottom).
#
# Same problem as ve_flow_fv.i (Nordbotten-Celia sloped aquifer) but the action defers
# geometry-variable declaration: the user declares z_top / z_bottom as MooseVariableFVReal
# and the action validates their existence + type, then couples to them. Mirrors the
# real-field workflow where geometry is read from the mesh. CSVDiff'd against the same
# nc_sloped_fv result, proving the opt-out path is equivalent.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 100
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 5.0e-5
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1020.0
    viscosity = 8.0e-4
  []
[]

[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    define_geometry_variables = false # user declares z_top / z_bottom below
  []
[]

# The user owns the geometry variables (as a mesh-read workflow would). They must be
# MooseVariableFVReal; the action checks this and errors otherwise.
[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * x'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x - 20'
  []
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []
  [nose_analytical]
    type = ParsedFunction
    expression = '1.0e-12 * 320.0 * 9.81 * 0.034899 * t / (0.2 * 5.0e-5)'
  []
[]

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_init
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.0
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

[FVBCs]
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [co2_injection]
    type = FVNeumannBC
    variable = sat_n
    boundary = left
    value = 2.8e-3
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
  end_time = 2.0e6
  dt = 1.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

[Postprocessors]
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_mass]
    type = FunctionValuePostprocessor
    function = '4.0e-6 * t / (0.2 * 20.0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [gravity_nose]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
