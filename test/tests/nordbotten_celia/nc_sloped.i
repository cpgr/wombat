# Nordbotten-Celia sloped-aquifer test -- Gamma >> 1 (gravity-dominated)
#
# Reference: Nordbotten & Celia (2012), Geological Storage of CO2, Chapter 3.
#
# Setup: 1D inclined aquifer, CO2 injected at x=0 (downdip end) via a
# DiracKernel point source. The formation dips updip in the +x direction:
#   z_top(x) = x * sin(2 deg) = 0.034899 * x  [m]
#
# The buoyancy term rho_c * |g| * grad(z_T) in the advective flux kernels
# drives CO2 updip (toward larger x). This test specifically exercises the
# sloped-surface gravity term that is zero in the flat BL test.
#
# Parameters:
#   K_up    = 1e-12 m2     phi = 0.20     H = 20 m
#   rho_n   = 700  kg/m3   mu_n = 5e-5 Pa*s
#   rho_w   = 1020 kg/m3   mu_w = 8e-4 Pa*s
#   theta   = 2 deg        sin(theta) = 0.034899
#   Q_inj   = 4e-6 m3/s   (volumetric CO2 injection rate)
#   Gamma   = K * (rho_w-rho_n) * g * sin(theta) * H / (Q * mu_n) = 11
#             (Gamma > 1: gravity-dominated regime)
#
# Analytical gravity nose velocity (Nordbotten & Celia, gravity-dominated limit):
#   v_N = K_up * (rho_w - rho_n) * g * sin(theta) / (phi * mu_n)
#       = 1e-12 * 320 * 9.81 * 0.034899 / (0.2 * 5e-5)
#       = 1.096e-5 m/s
#
# At t = 2e6 s: nose expected at x = v_N * t = 21.9 m.
#
# Verification metrics (postprocessors):
#   co2_integral:          integral(sat_n, 0, 100) -- mass conservation
#   analytical_mass:       Q_inj * t / (phi * H)  -- expected co2_integral
#   gravity_nose:          v_N * t                 -- expected nose position
#
# Mass conservation check: co2_integral should track analytical_mass closely.
# Nose position check: compare gravity_nose with LineValueSampler profile.

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

# ---------------------------------------------------------------------------
# Primary variables
# ---------------------------------------------------------------------------
[Variables]
  [pp_top]
  []
  [sat_n]
  []
[]

# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------
[Functions]
  # Hydrostatic brine pressure at top surface for a fully brine-saturated
  # formation. Sets pp_top = 0 at x = 100 m (open right boundary) and
  # increases linearly toward x = 0 (downdip, higher pressure).
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []

  # Linearly varying top elevation: formation dips updip in the +x direction.
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * x'
  []

  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x - 20'
  []

  # Analytical CO2 volume from mass conservation:
  #   integral(sat_n, 0, L) = Q_inj * t / (phi * H)
  [mass_analytical]
    type = ParsedFunction
    expression = '4.0e-6 * t / (0.2 * 20.0)'
  []

  # Analytical gravity nose position (gravity-dominated, Nordbotten-Celia):
  #   x_N(t) = K_up * (rho_w - rho_n) * g * sin(theta) / (phi * mu_n) * t
  [nose_analytical]
    type = ParsedFunction
    expression = '1.0e-12 * 320.0 * 9.81 * 0.034899 * t / (0.2 * 5.0e-5)'
  []
[]

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Geometry AuxVariables (tilted formation, held constant by no AuxKernel)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------
[Kernels]
  # CO2 mass equation (variable: sat_n, phase 0)
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

  # Brine mass equation (variable: pp_top, phase 1)
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

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------
[BCs]
  # Right boundary (updip end, x = 100 m): open outflow for brine.
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []

  # Left boundary (downdip, x = 0): CO2 injection as a constant mass flux.
  # PorousFlowSink with negative flux_function = source (injection into domain).
  # flux_function = -rho_n * Q_inj = -700 [kg/m3] * 4e-6 [m3/s] = -2.8e-3 [kg/m2/s].
  # use_mobility = use_relperm = multiply_by_density = false: the flux is taken
  # directly as specified, with no additional material-property multipliers.
  # This avoids requiring standard PorousFlow nodal properties that the VE
  # material system does not provide.
  # [co2_injection]
  #   type = PorousFlowSink
  #   PorousFlowDictator = ve_dictator
  #   variable = sat_n
  #   boundary = left
  #   fluid_phase = 0
  #   flux_function = -2.8e-3
  # []
  [co2_injection]
    type = NeumannBC
    variable = sat_n
    boundary = left
    value = 2.8e-3
  []
[]

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
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
  # Distinct densities: rho_w > rho_n drives updip CO2 migration.
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

# ---------------------------------------------------------------------------
# Executioner
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Postprocessors
# ---------------------------------------------------------------------------
[Postprocessors]
  # Mass conservation: integral of sat_n over [0, 100] m.
  # Should equal Q_inj * t / (phi * H) = 1e-6 * t [m].
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Analytical CO2 volume from injection mass balance.
  [analytical_mass]
    type = FunctionValuePostprocessor
    function = mass_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Analytical gravity nose position x_N = v_N * t.
  [gravity_nose]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Sample sat_n at monitoring points to track plume arrival.
  # Expected arrival times: t_arr(x) = x / v_N = x * phi * mu_n / (K * drho * g * sin(theta))
  #   x =  5 m -> t_arr = 4.56e5 s
  #   x = 10 m -> t_arr = 9.13e5 s
  #   x = 15 m -> t_arr = 1.37e6 s
  #   x = 20 m -> t_arr = 1.83e6 s
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

# ---------------------------------------------------------------------------
# VectorPostprocessors -- full saturation profile at end of each timestep
# ---------------------------------------------------------------------------
[VectorPostprocessors]
  [saturation_profile]
    type = LineValueSampler
    variable = sat_n
    start_point = '0 0 0'
    end_point = '100 0 0'
    num_points = 201
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
