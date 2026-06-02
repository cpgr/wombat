# Unit test for VEPlumeReconstruction in capillary_fringe mode with
# VECapillaryPressureTable.
#
# Flat formation: z_top = 0, z_bottom = -20 m, H = 20 m.
# Uniform sat_n = 0.1, uniform pp_top = 1e7 Pa.  Dirichlet BCs match ICs.
#
# Fluid properties: rho_co2 = 700, rho_brine = 1020 kg/m3
#   => delta_rho = 320 kg/m3
#   => Pc_max = delta_rho * g * H = 320 * 9.81 * 20 = 62784 Pa
#
# Linear Pc table: Sw(Pc) = 1 - 0.8 * Pc / 62784
#   pc_points = '0  15696  31392  47088  62784'
#   sw_points = '1.0  0.8  0.6  0.4  0.2'
#   (four equal intervals spanning [0, Pc_max])
#
# Analytical solution for h:
#   Sn(zeta) = 1 - Sw(Pc(zeta)) = 0.8 * delta_rho * g * zeta / Pc_max
#            = 0.8 * zeta / H  = 0.04 * zeta
#   sat_n_bar = (1/H) * integral_0^h Sn(zeta) dzeta
#             = (1/20) * 0.04 * h^2 / 2
#             = 0.001 * h^2
#   => h = sqrt(sat_n / 0.001) = sqrt(0.1 / 0.001) = sqrt(100) = 10.0 m
#
# Because Sn(zeta) is linear, the 64-point trapezoidal rule is exact,
# so Newton inversion returns h = 10.0 m to within convergence tolerance.

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

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 1.0e-3
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1020.0
    viscosity = 1.0e-3
  []
[]

[UserObjects]
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = capillary_fringe
  []
  [ve_pc_table]
    type = VECapillaryPressureTable
    sat_lr = 0.2
    pc_points = '0  15696  31392  47088  62784'
    sw_points = '1.0  0.8  0.6  0.4  0.2'
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
    value = 0.1
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
    initial_condition = -20.0
  []
  [h_plume]
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
    value = 0.1
  []
  [sat_right]
    type = DirichletBC
    variable = sat_n
    boundary = right
    value = 0.1
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
  [plume_reconstruction]
    type = VEPlumeReconstruction
    mode = capillary_fringe
    sat_n = sat_n
    S_wr = 0.2
    pc_uo = ve_pc_table
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

[Postprocessors]
  [h_avg]
    type = ElementAverageValue
    variable = h_plume
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
