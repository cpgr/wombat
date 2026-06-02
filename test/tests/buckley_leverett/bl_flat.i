# Buckley-Leverett verification test -- 1D flat formation, M = 1
#
# Equal density and viscosity for CO2 and brine (M = 1), linear relperms
# (S_wr = 0, krn_max = krw_max = 1), flat formation (z_top = 0, z_bottom = -10).
# No gravity contribution because grad(z_top) = 0 everywhere.
#
# With M = 1 the fractional flow function is f_n = S_n.  The saturation
# equation reduces to a linear advection PDE with wave speed:
#
#   v = K_up * |dp/dx| / (phi * mu)
#     = 1e-10 [m2] * 1000 [Pa/m] / (0.2 * 1e-3 [Pa*s])
#     = 5e-4 m/s
#
# Verification metric: integral of sat_n over [0, 100] m equals the front
# position x_front(t) = v * t for a step-function profile.
# At t = 1e5 s: x_front = 50 m, so the integral should be 50 m^2 / m = 50 m.
#
# Postprocessor 'co2_volume'    -- numerical  integral of sat_n
# Postprocessor 'analytical_front' -- 5e-4 * t  (front position)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
[]

# ---------------------------------------------------------------------------
# Pass VEDictator to all kernels and materials through GlobalParams.
# ---------------------------------------------------------------------------
[GlobalParams]
  VEDictator = ve_dictator
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
# Functions for ICs and postprocessors
# ---------------------------------------------------------------------------
[Functions]
  # Linear pressure from 1e5 Pa at x=0 down to 0 Pa at x=100 m.
  [pp_top_linear]
    type = ParsedFunction
    expression = '1.0e5 * (1.0 - x / 100.0)'
  []

  # Analytical BL front position as a function of time.
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

# ---------------------------------------------------------------------------
# Constant geometry AuxVariables (flat formation, no topography)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Kernels -- one storage + one flux kernel per primary variable
# ---------------------------------------------------------------------------
[Kernels]
  # CO2 (non-wetting, phase 0) mass equation -- acts on sat_n
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

  # Brine (wetting, phase 1) mass equation -- acts on pp_top
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
  # Pressure: Dirichlet on both ends to fix the pressure gradient.
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

  # Saturation: Dirichlet at the inlet (x = 0) only.
  # The right boundary is a natural outflow -- no BC required.
  [sat_n_left]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 1.0
  []
[]

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
[Materials]
  # Formation geometry (H = 10 m, flat top surface)
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
  []

  # Depth-averaged porosity
  [porosity]
    type = VEPorosity
    phi_bar = 0.2
  []

  # Upscaled horizontal permeability (isotropic 2x2, no H factor here)
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
  []

  # Constant fluid properties -- equal for both phases so M = 1
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
  []

  # Wraps sat_n primary variable into the ve_saturation material vector
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []

  # Sharp-interface upscaled relperms (linear: kr_n = S_n, kr_w = 1 - S_n)
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

# ---------------------------------------------------------------------------
# Postprocessors -- numerical vs analytical front position
# ---------------------------------------------------------------------------
[Postprocessors]
  # Integral of sat_n over the domain [0, 100] m.
  # For a step-function profile sat_n = 1 behind the front and 0 ahead,
  # this equals x_front numerically (units: m).
  [co2_volume]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Analytical front position x_front = v * t = 5e-4 * t [m].
  [analytical_front]
    type = FunctionValuePostprocessor
    function = front_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# ---------------------------------------------------------------------------
# VectorPostprocessors -- full saturation profile at every output step
# ---------------------------------------------------------------------------
[VectorPostprocessors]
  [saturation_profile]
    type = LineValueSampler
    variable = sat_n
    start_point = '0 0 0'
    end_point = '100 0 0'
    num_points = 101
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
