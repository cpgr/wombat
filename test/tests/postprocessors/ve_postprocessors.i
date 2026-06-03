# Unit test for all five VE diagnostic postprocessors.
#
# Same flat-formation setup as ve_auxkernels.i: uniform sat_n = 0.6,
# uniform pp_top = 1e7 Pa, H = 20 m, phi_bar = 0.2, rho_co2 = 700 kg/m3,
# rho_brine = 1000 kg/m3.  One timestep; fields remain at IC exactly.
#
# Expected postprocessor values:
#
#   footprint   = length where sat_n > 1e-6
#               = 1000 m  (full domain, since sat_n = 0.6 everywhere)
#
#   mass_mobile = rho_co2 * phi_bar * H * sat_n_mobile * domain_length
#               = 700 * 0.2 * 20 * 0.6 * 1000
#               = 1 680 000 kg
#
#   mass_trapped  = 0 kg  (sat_n_trap not coupled, defaults to zero)
#
#   mass_dissolved = 0 kg  (X_co2 not coupled, defaults to zero)
#
#   max_pp      = max(pp_top) = 1e7 Pa  (framework NodalMaxValue)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 1000
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
  end_time = 1.0
  dt = 1.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
[]

[Postprocessors]
  # VEPlumeFootprintPostprocessor:
  # sat_n = 0.6 > threshold = 1e-6 everywhere -> integral of 1 dA = 1000 m
  [footprint]
    type = VEPlumeFootprintPostprocessor
    sat_n = sat_n
    threshold = 1.0e-6
    execute_on = 'initial timestep_end'
  []

  # VEMobileCO2MassPostprocessor (no sat_n_trap coupled):
  # 700 * 0.2 * 20 * 0.6 * 1000 = 1 680 000 kg
  [mass_mobile]
    type = VEMobileCO2MassPostprocessor
    sat_n = sat_n
    rho_co2 = 700.0
    execute_on = 'initial timestep_end'
  []

  # VETrappedCO2MassPostprocessor (sat_n_trap defaults to 0):
  # 0 kg
  [mass_trapped]
    type = VETrappedCO2MassPostprocessor
    rho_co2 = 700.0
    execute_on = 'initial timestep_end'
  []

  # VEDissolvedCO2MassPostprocessor (X_co2 defaults to 0):
  # 0 kg
  [mass_dissolved]
    type = VEDissolvedCO2MassPostprocessor
    sat_n = sat_n
    rho_brine = 1000.0
    execute_on = 'initial timestep_end'
  []

  # NodalMaxValue on pp_top -- framework built-in, no custom class needed.
  # max(pp_top) = 1e7 Pa.  Subtract your reference pressure in post-processing
  # to obtain the overpressure for caprock-integrity screening.
  [max_pp]
    type = NodalExtremeValue
    variable = pp_top
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
