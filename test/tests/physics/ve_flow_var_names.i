# VE flow [Physics] -- renamed primary variables.
#
# Exercises the pressure_variable / saturation_variable rename parameters: instead of
# the default 'pp_top' / 'sat_n', the primary variables are named 'p_top' and 's_co2'.
# The Physics action must propagate the new names to every kernel and material it
# creates (variable couplings, kernel variable= parameters).
#
# Mirror of buckley_leverett/bl_flat.i (1D, M = 1, flat formation). The Physics block
# replaces the hand-written [Variables]/[Kernels]/[Materials], all user-facing objects
# (BCs, ICs, postprocessors) use the renamed variables, and the result is CSVDiff'd
# against the bl_flat gold -- proving the rename is transparent to the physics.

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
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Physics/VEFlow/FiniteElement]
  [ve]
    pressure_variable  = p_top   # renamed from default 'pp_top'
    saturation_variable = s_co2  # renamed from default 'sat_n'
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
  [p_top_ic]
    type = FunctionIC
    variable = p_top    # renamed
    function = pp_top_linear
  []
  [s_co2_ic]
    type = ConstantIC
    variable = s_co2    # renamed
    value = 0.0
  []
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
  [p_top_left]
    type = DirichletBC
    variable = p_top    # renamed
    boundary = left
    value = 1.0e5
  []
  [p_top_right]
    type = DirichletBC
    variable = p_top    # renamed
    boundary = right
    value = 0.0
  []
  [s_co2_left]
    type = DirichletBC
    variable = s_co2    # renamed
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
    variable = s_co2    # renamed; postprocessor name unchanged so gold matches
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
