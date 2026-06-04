# VE flow [Physics] -- FV, capillary = capillary_fringe, laterally varying H.
#
# Companion to ve_flow_fv_capillary_fringe.i (flat H): this one puts the FV fringe
# flux on a laterally varying thickness so the grad(H) fringe term
#   ve_dPcup_dH = delta_rho*|g|*sat_n/S_n(h)  (times grad(H).n)
# is ALSO exercised in the Jacobian, mirroring the FE ve_flow_fe_capillary_fringe.i
# coverage. This is the no-flow capillary-equilibrium geometry of
# capillary_equilibrium/cap_equilibrium_fv.i:
#   - flat top z_T = 0 (kills buoyancy), H(x) = 10 + 0.2 x.
#   - sat_n IC ~ a fringe no-flow equilibrium: with this LINEAR Sw(Pc) table the
#     column average is sat_n = 0.4*delta_rho*g*h^2/(Pc_max*H), so h=const gives
#     sat_n*H = const. sat_n = 1.2/H(x) is in [0.04, 0.12] (h ~ 8 m < H,
#     S_n(h) ~ 0.3) -- interior and off the plume-nose floor.
#
# Because this is a NO-FLOW state the Darcy velocity is ~0 everywhere, which sits
# exactly on the FV upwind-switch kink (HANDOVER Gotcha C). advected_interp_method
# = average (now a [Physics/VEFlow/FiniteVolume] option) selects the central,
# fully-differentiable interpolation so the capillary Jacobian can be checked here
# -- the capillary term itself is untouched by that choice.
#
# Verified by PetscJacobianTester (tests spec). Linear table => exact dSw/dPc =>
# strict 1e-7 ratio.

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
    type = VECapillaryPressureTable
    sat_lr = 0.2
    pc_points = '0  15696  31392  47088  62784'
    sw_points = '1.0  0.8  0.6  0.4  0.2'
  []
[]

[Physics/VEFlow/FiniteVolume]
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
    S_wr = 0.2
    advected_interp_method = average # differentiable at the no-flow (zero-velocity) state
    define_geometry_variables = false
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
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = '0.0'
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_fn
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
  []
[]

[FVBCs]
  # Single pressure pin removes the closed-domain null space (no sat_n BC: closed).
  [pp_pin]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 2.0e3
  dt = 2.0e3

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
  l_tol = 1.0e-8
[]
