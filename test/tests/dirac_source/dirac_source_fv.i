# Dirac point-source well -- FV formulation.
#
# Finite-volume counterpart of dirac_source.i: identical geometry, well, rate and
# checks, but with MooseVariableFVReal unknowns + VEFVMassTimeDerivative /
# VEFVAdvectiveFlux and the functor relperm. ConstantPointSource works on the FV sat_n
# variable too (adds the source to the cell containing the point).
#
# Same two checks as the FE version:
#   1. MASS CONSERVATION: co2_integral = integral(sat_n) tracks value*t/(rho_n*phi*H).
#   2. UPDIP MIGRATION: sat 20 m updip of the well rises while 20 m downdip stays ~0.
#
# FV upwinding is monotone, so unlike the FE version this carries the sharp plume
# without ringing. Geometry is set as FV variables by FunctionIC (live before the first
# Jacobian -- the FV idiom; see test/tests/cross_code or nc_sloped_fv).

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 20
  xmin = 0
  xmax = 200
  ymin = 0
  ymax = 100
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
    expression = '0.02 * x'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.02 * x - 20'
  []
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.02 * (200 - x)'
  []
  # Injected CO2 volume value*t/(rho_n*phi*H) plus the 1e-6 background (1e-6*200*100).
  [mass_analytical]
    type = ParsedFunction
    expression = '1.0e-6 * 200 * 100 + 0.3 * t / (700 * 0.2 * 20)'
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
  # Trace uniform background CO2 saturation. Injecting into a perfectly dry domain
  # (sat_n = 0) sits on the sharp-relperm kink where d(kr_n)/d(sat_n) clamps to 0
  # (VERelPermSharpUO uses x > 0 strictly; HANDOVER Gotcha A), dropping the CO2 mobility
  # from the Jacobian so Newton crawls (linear solve fine, nonlinear |R| barely moves).
  # Lifting the whole field just off zero restores the derivative. Adds 1e-6*area = 0.02
  # to co2_integral, folded into mass_analytical.
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 1.0e-6
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

[DiracKernels]
  # Interior CO2 well: 0.3 kg/s at (52.5, 52.5), on the FV sat_n (mass) equation.
  [co2_well]
    type = ConstantPointSource
    variable = sat_n
    point = '52.5 52.5 0'
    value = 0.3
  []
[]

[FVBCs]
  [pp_updip]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
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
  end_time = 3.0e6

  # Adaptive dt paces the sharp injection startup (sat jumps fast at the well as the
  # plume establishes), growing to dtmax once it settles. The CSV is synced to fixed
  # times (below), so the stepping does not affect the gold.
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e3
    growth_factor = 2.0
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 1.0e5

  solve_type = NEWTON
  # 'basic' = full Newton step (no line search). The default backtracking line search
  # ('bt') cuts the step to a tiny lambda during the sharp injection startup, stalling
  # convergence to a crawl (a perfect linear solve but a constant, tiny nonlinear |R|
  # decrement); the full step converges. Standard for PorousFlow-type problems.
  petsc_options_iname = '-pc_type -pc_factor_shift_type -snes_linesearch_type'
  petsc_options_value = 'lu NONZERO basic'

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
    function = mass_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_well]
    type = PointValue
    variable = sat_n
    point = '52.5 52.5 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_updip]
    type = PointValue
    variable = sat_n
    point = '72.5 52.5 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_downdip]
    type = PointValue
    variable = sat_n
    point = '32.5 52.5 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  # CSV gold lands on FIXED sync times, so it is deterministic (round times) regardless
  # of how the adaptive stepper subdivides -- the regression stays robust to convergence
  # changes. exodus still writes every step for visualisation.
  [csv]
    type = CSV
    sync_times = '3.0e5 6.0e5 9.0e5 1.2e6 1.5e6 1.8e6 2.1e6 2.4e6 2.7e6 3.0e6'
    sync_only = true
  []
[]
