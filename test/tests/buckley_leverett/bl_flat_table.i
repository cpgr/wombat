# Buckley-Leverett verification test -- tabulated relative permeability.
#
# Identical to bl_flat.i except the relative permeability is supplied by a
# tabulated VERelPermTableUO instead of the analytical VERelPermSharpUO. The
# table samples the SAME curve (kr_n = sat_n, kr_w = 1 - sat_n, i.e. S_wr = 0,
# krn_max = krw_max = 1) at several points on the line, so the piecewise-linear
# interpolant reproduces the analytical curve exactly. The solve therefore
# reproduces the bl_flat result (gold reused), verifying that the tabulated
# UserObject path reproduces analytical mode -- the field-case relperm plugs into
# the same FE/FV kernels with no other change.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  VEDictator = ve_dictator
[]

[UserObjects]
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = sharp_interface
  []
  # Tabulated relperm sampling kr_n = sat_n and kr_w = 1 - sat_n on the line.
  [relperm_uo]
    type = VERelPermTableUO
    sat_n_points = '0    0.25  0.5  0.75  1'
    krn_points   = '0    0.25  0.5  0.75  1'
    krw_points   = '1    0.75  0.5  0.25  0'
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
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
[]

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
    initial_condition = -10.0
  []
[]

[Kernels]
  [co2_storage]
    type = VEMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
  []
  [co2_flux]
    type = VEAdvectiveFluxS
    variable = sat_n
    fluid_phase = 0
    pp_top = pp_top
  []
  [brine_storage]
    type = VEMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
  []
  [brine_flux]
    type = VEAdvectiveFluxP
    variable = pp_top
    fluid_phase = 1
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

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
  []
  [porosity]
    type = VEPorosity
    phi_bar = 0.2
  []
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
  []
  [fluid_props]
    type = VEFluidPropertiesConst
    rho_co2 = 1000.0
    rho_brine = 1000.0
    mu_co2 = 1.0e-3
    mu_brine = 1.0e-3
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
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
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
