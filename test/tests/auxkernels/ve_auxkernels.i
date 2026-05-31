# Unit test for VEPlumeHeightAux and VEGravityNumberAux.
#
# A flat formation (z_top = 0, z_bottom = -20 m, H = 20 m) is initialised
# with a uniform depth-averaged CO2 saturation sat_n = 0.6 and uniform
# pressure pp_top = 1e7 Pa.  Dirichlet BCs match the ICs on both ends, so
# there is no pressure or saturation gradient -- the system is at steady
# state and the solver converges in a single Newton step.
#
# Expected AuxVariable values (verified via ElementAverageValue):
#
#   h_plume = sat_n * H / (1 - S_wr)
#           = 0.6 * 20 / (1 - 0.2)
#           = 15.0 m
#
#   gamma_ve = k_v * delta_rho * g * H^2 / (mu_n * phi_bar * Q * L)
#            = 1e-12 * 300 * 9.81 * 400 / (1e-3 * 0.2 * 1e-3 * 1000)
#            = 0.005886

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 1000
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
    S_wr = 0.2
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

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e7
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.6
  []
[]

# ---------------------------------------------------------------------------
# Geometry AuxVariables -- flat formation, no topography
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
    initial_condition = -20.0
  []
  [h_plume]
    order = CONSTANT
    family = MONOMIAL
  []
  [gamma_ve]
    order = CONSTANT
    family = MONOMIAL
  []
  [pc_up]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plume_height]
    type = VEPlumeHeightAux
    variable = h_plume
    execute_on = 'TIMESTEP_END'
  []
  [gravity_number]
    type = VEGravityNumberAux
    variable = gamma_ve
    k_v = 1.0e-12
    delta_rho = 300.0
    gravity = '0 0 -9.81'
    mu_n = 1.0e-3
    Q = 1.0e-3
    L = 1000.0
    execute_on = 'TIMESTEP_END'
  []
  [cap_pressure]
    type = ADMaterialRealAux
    variable = pc_up
    property = ve_pc_up
    execute_on = 'TIMESTEP_END'
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
    value = 1.0e7
  []
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 1.0e7
  []
  [sat_left]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 0.6
  []
  [sat_right]
    type = DirichletBC
    variable = sat_n
    boundary = right
    value = 0.6
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
    S_wr = 0.2
  []
  [cap_pressure]
    type = VEUpscaledCapPressure
    gravity = '0 0 -9.81'
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]

[Executioner]
  type = Transient
  end_time = 1.0
  dt = 1.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
[]

# ---------------------------------------------------------------------------
# Verify AuxVariable values via ElementAverageValue postprocessors.
# Both fields are uniform so the average equals the pointwise value.
# ---------------------------------------------------------------------------
[Postprocessors]
  [h_avg]
    type = ElementAverageValue
    variable = h_plume
    execute_on = 'TIMESTEP_END'
  []
  [gamma_avg]
    type = ElementAverageValue
    variable = gamma_ve
    execute_on = 'TIMESTEP_END'
  []
  # Pc^up = (rho_w - rho_n) * g * h = (1000 - 700) * 9.81 * 15 = 44145 Pa
  [pc_up_avg]
    type = ElementAverageValue
    variable = pc_up
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
