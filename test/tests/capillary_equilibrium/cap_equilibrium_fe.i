# Capillary no-flow equilibrium on a laterally varying thickness -- FE
#
# Purpose: verify the sat_n * grad(H) chain-rule part of the capillary drive
# grad(Pc^up) = d(Pc^up)/d(sat_n) * grad(sat_n) + d(Pc^up)/d(H) * grad(H).
#
# Setup (everything chosen so the ONLY active driver is the capillary term):
#   - flat top surface z_T = 0  -> the buoyancy drive rho*g*grad(z_T) is zero.
#   - laterally varying thickness H(x) = 10 + 0.2 x  (z_B = -H, so grad(z_T) = 0
#     but grad(H) = 0.2 /= 0).
#   - closed domain (no sat_n BC, single pp_top pin), no wells.
#   - constant per-phase densities; S_wr = 0.
#
# With the sharp-interface upscaled Pc^up = delta_rho * g * sat_n * H, a no-flow
# steady state requires grad(Pc^up) = 0, i.e. sat_n * H = const. The initial
# condition is set EXACTLY to that equilibrium:
#
#   sat_n_eq(x) = C / H(x),  C = 6   ->  sat_n in [0.2, 0.6], strictly interior.
#
# grad(Pc^up) = delta_rho * g * grad(sat_n * H) = delta_rho * g * grad(C) = 0,
# so the configuration is a fixed point: sat_n must NOT move.
#
# Discriminating power: WITHOUT the grad(H) term the kernel would instead drive
# flux ~ d(Pc^up)/d(sat_n) * grad(sat_n) = delta_rho*g*H*grad(C/H) /= 0, pushing
# sat_n toward uniform and away from the prescribed profile. The ElementL2Error
# satn_drift therefore stays ~O(1e-4) WITH the term and grows to O(1e-2)+ without
# it (which also flips the would-be steady state). This checks the SIGN and
# MAGNITUDE of the grad(H) term, which the Jacobian test cannot.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
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

[Functions]
  [z_bottom_fn]
    type = ParsedFunction
    expression = '-(10.0 + 0.2 * x)'
  []
  [sat_n_eq]
    type = ParsedFunction
    expression = '6.0 / (10.0 + 0.2 * x)'
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e5
  []
  [sat_n_ic]
    type = FunctionIC
    variable = sat_n
    function = sat_n_eq
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_fn
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
  # Single pressure pin removes the closed-domain pressure null space. Its value
  # equals the uniform IC, so the equilibrium pp_top = 1e5 satisfies it exactly.
  [pp_pin]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
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
  end_time = 2.0e4
  dt = 2.0e3

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  [satn_drift]
    type = ElementL2Error
    variable = sat_n
    function = sat_n_eq
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
