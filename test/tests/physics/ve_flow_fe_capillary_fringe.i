# VE flow [Physics] -- FE, capillary = capillary_fringe.
#
# Physics-generated fringe counterpart of ve_flow_fe_capillary.i. capillary =
# capillary_fringe makes the FE physics run VEPlumeReconstruction in Newton-
# inversion mode against the shared pc_uo Sw(Pc) table, and makes
# VEUpscaledCapPressure use the fringe coefficients
#   ve_dPcup_dsatn = delta_rho*|g|*H/S_n(h),  ve_dPcup_dH = delta_rho*|g|*sat_n/S_n(h)
# (vs the sharp H/(1-S_wr), sat_n/(1-S_wr)), with S_n(h) = 1 - Sw(delta_rho*|g|*h).
#
# Geometry/state (mirrors capillary_equilibrium/cap_equilibrium_fe.i so BOTH the
# grad(sat_n) and grad(H) fringe terms are active at a NON-uniform interior state,
# which the AD-consistent fringe Jacobian must match):
#   - flat top z_T = 0 (kills buoyancy), laterally varying H(x) = 10 + 0.2 x.
#   - sat_n IC ~ a fringe no-flow equilibrium: with this LINEAR Sw(Pc) table the
#     column average is sat_n = 0.4*delta_rho*g*h^2/(Pc_max*H), so h=const gives
#     sat_n*H = const, the same form as the sharp case. sat_n = 1.2/H(x) lies in
#     [0.04, 0.12] (thin plume, h ~ 8 m < H, S_n(h) ~ 0.3) -- strictly interior and
#     away from the plume-nose floor, so no relperm / S_n kink (cf. Gotchas A, D).
#
# Verified by PetscJacobianTester (tests spec). The linear table makes
# dSw/dPc piecewise-constant (exact), so the strict 1e-7 ratio holds (Gotcha-D
# analog: a nonlinear table would amplify the value-vs-derivative inconsistency).

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
    density = 1000.0 # delta_rho = 300 kg/m3
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
  [ve_pc_table]
    # Linear Sw(Pc): Sw = 1 - 0.8 * Pc / 62784, sat_lr = 0.2. dSw/dPc is exact
    # (piecewise constant), so the fringe Jacobian matches FD to 1e-7.
    type = VECapillaryPressureTable
    sat_lr = 0.2
    pc_points = '0  15696  31392  47088  62784'
    sw_points = '1.0  0.8  0.6  0.4  0.2'
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
    capillary = capillary_fringe
    pc_uo = ve_pc_table
    S_wr = 0.2 # matches the table sat_lr (used for the Newton warm start)
    define_geometry_variables = false # z_top/z_bottom declared in [AuxVariables] below
  []
[]

[Functions]
  [z_bottom_fn]
    type = ParsedFunction
    expression = '-(10.0 + 0.2 * x)'
  []
  [sat_n_eq]
    type = ParsedFunction
    expression = '1.2 / (10.0 + 0.2 * x)'
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

[BCs]
  # Single pressure pin removes the closed-domain null space (no sat_n BC: closed).
  [pp_pin]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
  []
[]

[Executioner]
  type = Transient
  end_time = 2.0e3
  dt = 2.0e3

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
  l_tol = 1.0e-8
[]
