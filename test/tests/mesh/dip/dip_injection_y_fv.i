# FV version of dip_injection_y.i -- y-dip only (90-deg rotation of x case).
#
# Finite-volume counterpart: MooseVariableFVReal unknowns, VEFVMassTimeDerivative
# + VEFVAdvectiveFlux, functor relperm. z_top = 0.034899*y (y-dip only); CO2
# injected through the centre 10 m of the bottom boundary (20 <= x <= 30) migrates
# updip in +y. co2_integral = value*L_inj*t/(H*phi*rho_n) = 1e-5*t (L_inj = 10 m).

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
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - y)'
  []
  [nose_analytical]
    type = ParsedFunction
    expression = '1.0e-12 * 320.0 * 9.81 * 0.034899 * t / (0.2 * 5.0e-5)'
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
  # Open top boundary (updip end, y = 100): Dirichlet brine pressure.
  [pp_top_bc]
    type = FVDirichletBC
    variable = pp_top
    boundary = top
    value = 0.0
  []
  # CO2 injection through the centre 10 m of the bottom (downdip) boundary.
  [co2_injection]
    type = FVNeumannBC
    variable = sat_n
    boundary = injection_zone
    value = 2.8e-3
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
[]

[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n = sat_n
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
  end_time = 5.0e6
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
  [gravity_nose_y]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_25m]
    type = PointValue
    variable = sat_n
    point = '25 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_50m]
    type = PointValue
    variable = sat_n
    point = '25 50 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mass conservation: value * L_inj * t / (H * phi * rho_n), L_inj = 10 m.
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
