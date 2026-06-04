# VE flow [Physics] -- FV, capillary = capillary_fringe.
#
# Physics-generated fringe counterpart of ve_flow_fv_capillary.i. capillary =
# capillary_fringe makes the FV physics build VEFVCapPressure in Newton-inversion
# mode against the shared pc_uo Sw(Pc) table, so the flux coefficients use the
# fringe S_n(h) = 1 - Sw(delta_rho*|g|*h) instead of the sharp (1 - S_wr).
#
# State: flat formation H = 10 m, a pressure gradient drives a NON-zero Darcy
# velocity (so the upwind branch is differentiable -- cf. Gotcha C, which bites
# the FV Jacobian only at a no-flow / zero-velocity state), and a uniform interior
# sat_n = 0.1 that the LINEAR Sw(Pc) table can represent on a 10 m column
# (max column-average ~0.19 here). This is the proven FV capillary-Jacobian
# pattern (ve_flow_fv_capillary.i) with the fringe table swapped in; it exercises
# the fringe Newton inversion + S_n AD propagation. The grad(H) term vanishes for
# flat H, as on the sharp FV Jacobian test.
#
# Verified by PetscJacobianTester (tests spec). Linear table => exact dSw/dPc =>
# strict 1e-7 ratio holds (Gotcha-D analog).

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
    type = VECapillaryPressureTableUO
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
    S_wr = 0.2 # matches the table sat_lr (used for the Newton warm start)
    define_geometry_variables = false # z_top/z_bottom declared in [AuxVariables] below
  []
[]

[Functions]
  [pp_top_linear]
    type = ParsedFunction
    expression = '1.0e5 * (1.0 - x / 100.0)'
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
    value = 0.1 # interior, representable by the table on a 10 m column
  []
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = '0.0'
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = '-10.0'
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
  [pp_left]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
  []
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [sat_n_left]
    type = FVDirichletBC
    variable = sat_n
    boundary = left
    value = 0.1
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
  end_time = 1.0e3
  dt = 1.0e3

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
  l_tol = 1.0e-8
[]
