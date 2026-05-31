# Nordbotten-Celia sloped-aquifer test -- y-dip only (90-degree rotation of
# dip_injection_x.i).
#
# Geometry:
#   z_top    =  0.034899 * y   (sin 2 deg in y only; no x-dip)
#   z_bottom = z_top - 20 m   (H = 20 m, uniform)
#   Domain   :  50 m x 100 m  (axes swapped relative to x-dip test)
#
# Injection: NeumannBC at the centre 10 m of the bottom (downdip) boundary:
# 20 m <= x <= 30 m at y = 0.  CO2 migrates updip in the +y direction only.
# With no x-dip the plume is symmetric about x = 25 m.
#
# Analytical y-nose velocity (Nordbotten-Celia, gravity-dominated):
#   v_Ny = K * (rho_w - rho_n) * g * sin_y / (phi * mu_n)
#        = 1e-12 * 320 * 9.81 * 0.034899 / (0.2 * 5e-5)
#        = 1.096e-5 m/s     (identical to dip_injection_x.i by symmetry)
#
# At t = 5e6 s: nose expected near y = 54.8 m (past the halfway mark of 50 m).

[Mesh]
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 100
    xmin = 0
    xmax = 50
    ymin = 0
    ymax = 100
  []
  # Injection zone: centre 10 m of the bottom boundary (20 <= x <= 30).
  [injection_zone]
    type = ParsedGenerateSideset
    input = base_mesh
    combinatorial_geometry = 'y < 1e-6 & x >= 20.0 & x <= 30.0'
    normal = '0 -1 0'
    new_sideset_name = 'injection_zone'
  []
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
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * y'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * y - 20'
  []
  # Hydrostatic brine pressure: pp_top = 0 at the open top boundary (y = 100).
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - y)'
  []
  # Analytical y-nose position (gravity-dominated Nordbotten-Celia).
  [nose_analytical]
    type = ParsedFunction
    expression = '1.0e-12 * 320.0 * 9.81 * 0.034899 * t / (0.2 * 5.0e-5)'
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
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

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
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
  # Open top boundary (updip end, y = 100): brine exits freely.
  [pp_top_bc]
    type = DirichletBC
    variable = pp_top
    boundary = top
    value = 0.0
  []
  [sat_top]
    type = DirichletBC
    variable = sat_n
    boundary = top
    value = 0.0
  []
  # CO2 injection through the centre 10 m of the bottom (downdip) boundary.
  # rho_n * Q_inj = 700 * 4e-6 = 2.8e-3 kg/(m*s).
  [co2_injection]
    type = NeumannBC
    variable = sat_n
    boundary = injection_zone
    value = 2.8e-3
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
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
  []
  [fluid_props]
    type = VEFluidPropertiesConst
    rho_co2 = 700.0
    rho_brine = 1020.0
    mu_co2 = 5.0e-5
    mu_brine = 8.0e-4
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
  end_time = 5.0e6
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
  # Analytical y-nose position.
  [gravity_nose_y]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Plume centreline: x = 25 m, y = 25 m (quarter-way).
  [sat_at_25m]
    type = PointValue
    variable = sat_n
    point = '25 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Plume centreline: x = 25 m, y = 50 m (halfway).
  [sat_at_50m]
    type = PointValue
    variable = sat_n
    point = '25 50 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mass conservation: the domain integral of sat_n equals the injected volume
  # value * L_inj * t / (H * phi * rho_n), with L_inj = 10 m (20 <= x <= 30).
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_mass]
    type = FunctionValuePostprocessor
    function = '2.8e-3 * 10.0 * t / (20.0 * 0.2 * 700.0)'
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
