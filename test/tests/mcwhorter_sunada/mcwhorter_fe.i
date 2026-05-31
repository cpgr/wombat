# McWhorter-Sunada capillary verification -- 1D counter-current imbibition, FE
#
# Validates the SIGN and MAGNITUDE of the VE capillary-diffusion term (which the
# Jacobian tests cannot). With per-phase-constant densities, a flat formation
# (grad z_T = 0) and a CLOSED far (right) boundary, the total volume flux
# q_t = q_n + q_w is uniform in x and pinned to zero -> pure counter-current
# capillary imbibition. The CO2 saturation then obeys the nonlinear diffusion
# equation
#
#   phi dS/dt = d/dx [ D(S) dS/dx ],   D(S) = (K * C / mu) * S * (1 - S),
#   C = (rho_w - rho_n) * g * H / (1 - S_wr) = d Pc^up / dS.
#
# The self-similar (Boltzmann) solution S(x/sqrt(t)) is produced by
# mcwhorter_reference.py -> mcwhorter_reference.csv. The ElementL2Error
# postprocessor co2_l2_error measures || sat_n - S_ref || at t = 1e5 s.
#
# A wrong capillary sign makes D < 0 (anti-diffusion): the solve diverges or the
# error explodes, so this test genuinely checks the sign.
#
# Parameters MUST stay in sync with mcwhorter_reference.py.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
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
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
  []
[]

[Functions]
  # Semi-analytical McWhorter-Sunada profile at t = 1e5 s (x, S_ref).
  # axis = x makes this a function of the x-coordinate (not time, the default).
  [mcwhorter_ref]
    type = PiecewiseLinear
    data_file = mcwhorter_reference.csv
    format = columns
    axis = x
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e6
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.05
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
    capillary = true
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

# Counter-current imbibition BCs:
#   left  (x=0): fixed CO2 saturation inlet + fixed brine pressure (lets brine
#                escape as CO2 imbibes -> q_w(0) = -q_n(0)).
#   right (x=L): no BC on either variable -> natural no-flow -> closes the
#                domain so q_t = 0 everywhere.
[BCs]
  [sat_n_inlet]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 0.7
  []
  [pp_top_left]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
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
    rho_co2 = 700.0
    rho_brine = 1000.0
    mu_co2 = 1.0e-3
    mu_brine = 1.0e-3
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [plume_reconstruction]
    type = VEPlumeReconstruction
    sat_n = sat_n
    S_wr = 0.0
  []
  [cap_pressure]
    type = VEUpscaledCapPressure
    gravity = '0 0 -9.81'
    S_wr = 0.0
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
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
  end_time = 1.0e5
  dt = 500.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-9
  nl_abs_tol = 1.0e-11
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  # L2 norm of (sat_n - analytical) at the final time. Near zero = match.
  [co2_l2_error]
    type = ElementL2Error
    variable = sat_n
    function = mcwhorter_ref
    execute_on = 'FINAL'
  []
[]

[VectorPostprocessors]
  # Full profile for plotting / debugging (not golded).
  [saturation_profile]
    type = LineValueSampler
    variable = sat_n
    start_point = '0 0 0'
    end_point = '100 0 0'
    num_points = 201
    sort_by = x
    execute_on = 'FINAL'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
