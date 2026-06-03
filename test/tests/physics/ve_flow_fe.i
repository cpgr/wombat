# VE flow [Physics] regression -- FE (continuous Galerkin).
#
# Reproduces nc_sloped.i (Nordbotten-Celia sloped aquifer) but replaces the
# hand-written [Variables], [Kernels] and [Materials] blocks with a single
# [Physics/VE/Flow/FiniteElement] block. The CSV postprocessor output is
# CSVDiff'd against the nc_sloped gold to prove the generated objects are
# numerically identical to the hand-assembled ones.
#
# Everything outside the flow setup -- mesh, fluid-property UOs, relperm UO,
# functions, ICs, the z_top/z_bottom geometry aux variables, the injection /
# outflow BCs, executioner and postprocessors -- stays the user's responsibility.

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

# All of [Variables], [Kernels] and [Materials] are generated here.
[Physics/VEFlow/FiniteElement]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
  []
[]

[Functions]
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * x'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x - 20'
  []
  [mass_analytical]
    type = ParsedFunction
    expression = '4.0e-6 * t / (0.2 * 20.0)'
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

# z_top / z_bottom are declared by the [Physics] action (FE LAGRANGE); the ICs above
# set their values. No [AuxVariables] block is needed for the geometry.

[BCs]
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [co2_injection]
    type = NeumannBC
    variable = sat_n
    boundary = left
    value = 2.8e-3
  []
[]

[Executioner]
  type = Transient
  end_time = 2.0e6
  dt = 1.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

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
    function = mass_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [gravity_nose]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_5m]
    type = PointValue
    variable = sat_n
    point = '5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_10m]
    type = PointValue
    variable = sat_n
    point = '10 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_20m]
    type = PointValue
    variable = sat_n
    point = '20 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
