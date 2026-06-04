# VE flow [Physics] -- FE, capillary = sharp_interface.
#
# Physics-generated counterpart of bl_flat_capillary_fe.i. The single
# capillary = sharp_interface flag makes the FE physics set capillary on VEAdvectiveFluxS
# AND create the VEPlumeReconstruction + VEUpscaledCapPressure materials that
# supply ve_dPcup_dsatn / ve_dPcup_dH. Verified by PetscJacobianTester (tests
# spec): the AD Jacobian of the generated capillary-diffusion term must match FD.

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
    capillary = sharp_interface
    S_wr = 0.0
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
