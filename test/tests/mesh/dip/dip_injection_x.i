# Nordbotten-Celia sloped-aquifer test -- x-dip only, injection at centre of
# left boundary.
#
# Geometry:
#   z_top    =  0.034899 * x   (sin 2 deg in x only; no y-dip)
#   z_bottom = z_top - 20 m   (H = 20 m, uniform)
#   Domain   :  100 m x 50 m
#
# Injection: NeumannBC at the middle 10 m section of the left (downdip)
# boundary: 20 m <= y <= 30 m.  CO2 migrates updip in the +x direction only.
# With no y-dip the plume is symmetric about y = 25 m.
#
# Analytical x-nose velocity (Nordbotten-Celia, gravity-dominated):
#   v_Nx = K * (rho_w - rho_n) * g * sin_x / (phi * mu_n)
#        = 1e-12 * 320 * 9.81 * 0.034899 / (0.2 * 5e-5)
#        = 1.096e-5 m/s
#
# At t = 5e6 s: nose expected near x = 54.8 m (past the halfway mark of 50 m).

[Mesh]
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 50
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 50
  []
  # Injection zone: centre 10 m of the left boundary (20 <= y <= 30).
  [injection_zone]
    type = ParsedGenerateSideset
    input = base_mesh
    combinatorial_geometry = 'x < 1e-6 & y >= 20.0 & y <= 30.0'
    normal = '-1 0 0'
    new_sideset_name = 'injection_zone'
  []
[]

[GlobalParams]
  VEDictator = ve_dictator
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
    expression = '0.034899 * x'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x - 20'
  []
  # Hydrostatic brine pressure: pp_top = 0 at the open right boundary (x = 100).
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []
  # Analytical x-nose position (gravity-dominated Nordbotten-Celia).
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
  # Open right boundary (updip end): brine exits freely.
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [sat_right]
    type = DirichletBC
    variable = sat_n
    boundary = right
    value = 0.0
  []
  # CO2 injection through the centre 10 m of the left (downdip) boundary.
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
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
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
  # Analytical x-nose position.
  [gravity_nose_x]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Plume centreline: x = 25 m (quarter-way), y = 25 m.
  [sat_at_25m]
    type = PointValue
    variable = sat_n
    point = '25 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Plume centreline: x = 50 m (halfway), y = 25 m.
  [sat_at_50m]
    type = PointValue
    variable = sat_n
    point = '50 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mass conservation: the domain integral of sat_n equals the injected volume
  # value * L_inj * t / (H * phi * rho_n), with L_inj = 10 m (20 <= y <= 30).
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
