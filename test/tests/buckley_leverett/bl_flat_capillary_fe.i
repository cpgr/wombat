# Buckley-Leverett with capillary pressure -- 1D flat formation, FE (two-pressure VE)
#
# Same geometry and fluid mobility as bl_flat.i but with:
#   - rho_co2 = 700 kg/m3, rho_brine = 1000 kg/m3  (delta_rho = 300 kg/m3)
#   - capillary = true on VEAdvectiveFluxS
#   - VEPlumeReconstruction + VEUpscaledCapPressure in [Materials]
#
# With linear relperms (M=1) and the sharp-interface capillary pressure
#   Pc^up = delta_rho * g * sat_n * H / (1 - S_wr)
#         = 300 * 9.81 * sat_n * 10 / 1   [S_wr=0]
#         = 29430 * sat_n  Pa
# the CO2 mass equation gains the capillary-diffusion term
#   grad(Pc^up) = 29430 * grad(sat_n)
# which broadens the plume front relative to the capillary-off case.
#
# Verification:
#   1. Mass is still conserved (co2_volume = analytical_front), though this
#      alone does NOT validate the capillary sign -- see McWhorter-Sunada.
#   2. PetscJacobianTester (see tests spec) validates the AD Jacobian.

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

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
