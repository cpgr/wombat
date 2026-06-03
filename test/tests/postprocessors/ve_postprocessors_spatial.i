# Postprocessor test with spatial saturation variation, non-zero trapping,
# and non-zero dissolution.  The plume occupies only part of the domain so
# VEPlumeFootprintPostprocessor returns a value less than the domain length.
#
# sat_n(x) = max(0, 0.6 - 0.001*x)
#   -> 0.6 at x=0, linear decrease, zero for x >= 600 m.
#
# sat_n_trap = 0.10  (constant AuxVariable, simulates VEHysteresis output)
# X_co2      = 0.05  (constant AuxVariable, simulates VEDissolution output)
#
# Only VEMassTimeDerivative kernels are active (no advective flux), so
# sat_n stays at its initial condition.  The initial residual is identically
# zero and Newton converges without a linear solve.
#
# Mesh: nx=10, xmin=0, xmax=1000 m -- elements [0,100], [100,200], ...
# H = 20 m, phi_bar = 0.2, rho_co2 = 700 kg/m3, rho_brine = 1000 kg/m3
#
# Q1 elements integrate linear functions exactly (2-point Gauss).
#
# Analytical expected values:
#
#   integral(sat_n, 0, 1000)
#     = (0.6 + 0)/2 * 600  +  0 * 400  =  180 m
#
#   integral(max(0, sat_n - sat_n_trap), 0, 1000)
#     clamp activates for x > 500 m (where sat_n < 0.1 = sat_n_trap)
#     = integral(0.5 - 0.001*x, 0, 500)
#     = [0.5x - x^2/2000] from 0 to 500
#     = 250 - 125 = 125 m
#
#   integral(1 - sat_n, 0, 1000)
#     = integral(0.4 + 0.001*x, 0, 600)  +  integral(1, 600, 1000)
#     = [0.4x + x^2/2000] from 0 to 600  +  400
#     = (240 + 180)  +  400  =  820 m
#
#   footprint    = 600 m
#     All elements [0,600) have sat_n > 1e-6 at all Gauss points.
#     Elements [600,1000] have sat_n = 0 everywhere.
#
#   mass_mobile  = 700 * 0.2 * 20 * 125  =  350 000 kg
#
#   mass_trapped = 700 * 0.2 * 20 * 0.1 * 1000  =  280 000 kg
#
#   mass_dissolved = 1000 * 0.2 * 20 * 0.05 * 820  =  164 000 kg
#
#   max_pp       = 1e7 Pa

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

[Functions]
  [sat_n_func]
    type = ParsedFunction
    expression = 'max(0.0, 0.6 - 0.001 * x)'
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
    type = FunctionIC
    variable = sat_n
    function = sat_n_func
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
  # Held at 0.10 throughout -- simulates the output of a future VEHysteresis.
  [sat_n_trap]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
  []
  # Held at 0.05 throughout -- simulates the output of a future VEDissolution.
  [X_co2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  []
[]

# Mass storage only -- no advective flux.  Initial residual is zero so
# Newton converges without a linear solve and fields stay at their ICs.
[Kernels]
  [co2_storage]
    type = VEMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
  []
  [brine_storage]
    type = VEMassTimeDerivative
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
[]

[Executioner]
  type = Transient
  end_time = 1.0
  dt = 1.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_abs_tol = 1.0e-12
  nl_rel_tol = 1.0e-10
[]

[Postprocessors]
  # 6 elements (0-600 m) have sat_n > 1e-6; 4 elements (600-1000 m) have
  # sat_n = 0. footprint = 6 * 100 m = 600 m.
  [footprint]
    type = VEPlumeFootprintPostprocessor
    sat_n = sat_n
    threshold = 1.0e-6
    execute_on = 'initial timestep_end'
  []

  # sat_n_mobile = max(0, sat_n - 0.1).  Clamp activates for x > 500 m.
  # mass_mobile = 700 * 0.2 * 20 * 125 = 350 000 kg.
  [mass_mobile]
    type = VEMobileCO2MassPostprocessor
    sat_n = sat_n
    sat_n_trap = sat_n_trap
    rho_co2 = 700.0
    execute_on = 'initial timestep_end'
  []

  # mass_trapped = 700 * 0.2 * 20 * 0.1 * 1000 = 280 000 kg.
  [mass_trapped]
    type = VETrappedCO2MassPostprocessor
    sat_n_trap = sat_n_trap
    rho_co2 = 700.0
    execute_on = 'initial timestep_end'
  []

  # mass_dissolved = 1000 * 0.2 * 20 * 0.05 * 820 = 164 000 kg.
  [mass_dissolved]
    type = VEDissolvedCO2MassPostprocessor
    sat_n = sat_n
    X_co2 = X_co2
    rho_brine = 1000.0
    execute_on = 'initial timestep_end'
  []

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
