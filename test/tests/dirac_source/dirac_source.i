# Dirac point-source well -- FE formulation.
#
# CO2 injected at an INTERIOR point (52.5, 52.5) via a ConstantPointSource DiracKernel
# (value in kg/s, since the VE mass equation storage H*phi*rho*S is in mass units; a
# positive value is a source). This is the first test of a point/well source -- every
# other injection in the suite is a boundary Neumann flux.
#
# The aquifer dips updip in +x (z_top = 0.02*x), so buoyancy drives the injected CO2
# updip. Two checks:
#   1. MASS CONSERVATION: co2_integral = integral(sat_n) must track the analytic
#      injected volume value*t/(rho_n*phi*H) (no CO2 leaves -- sat_n has no outflow BC
#      and the plume stays interior). This verifies the Dirac source rate/units.
#   2. UPDIP MIGRATION: sat at a point 20 m UPDIP of the well rises while a point 20 m
#      DOWNDIP stays ~0 -- the buoyant tongue extends updip, not downdip.
#
# FE caveat: continuous-Galerkin advection of the sharp point-source plume can ring
# (HANDOVER Phase 4b note 2); dirac_source_fv.i is the robust counterpart. Parameters
# mirror nordbotten_celia/nc_sloped (K=1e-12, phi=0.2, H=20).

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
  # Trace uniform background CO2 saturation. Injecting into a perfectly dry domain
  # (sat_n = 0) sits on the sharp-relperm kink where d(kr_n)/d(sat_n) clamps to 0
  # (VERelPermSharpUO uses x > 0 strictly; HANDOVER Gotcha A), dropping the CO2 mobility
  # from the Jacobian so Newton crawls. Lifting the whole field just off zero restores
  # the derivative. Adds 1e-6*area = 0.02 to co2_integral, folded into mass_analytical.
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

[DiracKernels]
  # Interior CO2 well: 0.3 kg/s injected at (52.5, 52.5). ConstantPointSource adds
  # -value*test, so a positive value is a source on the sat_n (mass) equation.
  [co2_well]
    type = ConstantPointSource
    variable = sat_n
    point = '52.5 52.5 0'
    value = 0.3
  []
[]

[BCs]
  # Open updip boundary: brine pressure Dirichlet. CO2 (sat_n) is natural no-flow; the
  # plume stays interior, so CO2 is conserved (mass-balance check).
  [pp_updip]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
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

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

[Postprocessors]
  # Mass conservation: co2_integral should track mass_analytical (the injected volume).
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
  # Updip migration: sat 20 m updip rises; 20 m downdip stays ~0.
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
