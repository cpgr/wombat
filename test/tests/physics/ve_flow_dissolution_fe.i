# VE flow [Physics] -- FE, enable_dissolution = true.
#
# Physics-generated counterpart of dissolution/constant_rate.i. enable_dissolution =
# true creates the FULL dissolution chain -- the VEDissolution material, the c_diss
# aux variable, the VEDissolvedCO2Aux accumulator, and the VEDissolutionSink kernel --
# from the dissolution_flux / dissolution_s_ref parameters. The geometry aux variables
# (z_top, z_bottom) are also action-declared; only the ICs/BCs remain user-written.
# CSVDiff against the constant_rate result.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 1.0e-3
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1000.0
    viscosity = 1.0e-3
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

[Physics/VEFlow/FiniteElement]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    enable_dissolution = true
    dissolution_flux = 1.4
    dissolution_s_ref = 0.05
  []
[]

# Only ICs/BCs are user-written; the action declares pp_top, sat_n, z_top, z_bottom and
# c_diss. The geometry elevations are given constant values here (flat 10 m formation).
[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e6
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5
  []
  [z_top_ic]
    type = ConstantIC
    variable = z_top
    value = 0.0
  []
  [z_bottom_ic]
    type = ConstantIC
    variable = z_bottom
    value = -10.0
  []
[]

[BCs]
  [pp_left]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
  []
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 1.0e6
  []
[]

[Executioner]
  type = Transient
  end_time = 100.0
  dt = 20.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
[]

[Postprocessors]
  [sat_n_avg]
    type = ElementAverageValue
    variable = sat_n
    execute_on = 'timestep_end final'
  []
  [dissolved_co2_mass]
    type = ElementIntegralVariablePostprocessor
    variable = c_diss
    execute_on = 'timestep_end final'
  []
  [mobile_co2_mass]
    type = VEMobileCO2MassPostprocessor
    sat_n = sat_n
    rho_co2 = 700.0
    execute_on = 'timestep_end final'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
