# VE flow [Physics] regression -- tabulated relative permeability (FE).
#
# Mirrors buckley_leverett/bl_flat_table.i, but replaces the hand-written
# [Variables]/[Kernels]/[Materials] blocks with a single
# [Physics/VEFlow/FiniteElement] block. The relperm is supplied by a tabulated
# VERelPermTableUO whose table samples the analytical sharp-interface curve
# (kr_n = sat_n, kr_w = 1 - sat_n, i.e. S_wr = 0, krn_max = krw_max = 1) at
# points on the line, so the piecewise-linear interpolant reproduces the
# analytical curve exactly. The solve therefore reproduces the bl_flat
# Buckley-Leverett result, proving the tabulated UO plugs into the Physics-
# generated FE kernels with no other change (the field-case relperm path).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 1000.0
    viscosity = 1.0e-3
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1000.0
    viscosity = 1.0e-3
  []
[]

[UserObjects]
  # Tabulated relperm sampling kr_n = sat_n and kr_w = 1 - sat_n on the line.
  [relperm_uo]
    type = VERelPermTableUO
    sat_n_points = '0    0.25  0.5  0.75  1'
    krn_points   = '0    0.25  0.5  0.75  1'
    krw_points   = '1    0.75  0.5  0.25  0'
  []
[]

# All of [Variables], [Kernels] and [Materials] are generated here.
[Physics/VEFlow/FiniteElement]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
  []
[]

[Functions]
  [pp_top_linear]
    type = ParsedFunction
    expression = '1.0e5 * (1.0 - x / 100.0)'
  []
  [front_analytical]
    type = ParsedFunction
    expression = '5.0e-4 * t'
  []
[]

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_linear
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.0
  []
  # z_top / z_bottom are declared by the [Physics] action (FE LAGRANGE); set their
  # constant values here. No [AuxVariables] block is needed for the geometry.
  [z_top_ic]
    type = ConstantIC
    variable = z_top
    value = 0.0
  []
  [z_bottom_ic]
    type = ConstantIC
    variable = z_bottom
    value = -10.0
  []
[]

[BCs]
  [pp_left]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
  []
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [sat_n_left]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 1.0
  []
[]

[Executioner]
  type = Transient
  end_time = 1.0e5
  dt = 1000.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  [co2_volume]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_front]
    type = FunctionValuePostprocessor
    function = front_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
