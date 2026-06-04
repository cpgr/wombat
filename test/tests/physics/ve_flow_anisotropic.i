# VE flow [Physics] -- FE, anisotropic K_up (off-diagonal K_up_xy).
#
# Physics-generated counterpart of permeability/anisotropic_permeability.i: a
# method-of-manufactured-solution check that K_up_xy is forwarded into the
# VEPermeability tensor and used by the Darcy operator. With
# K_up = [[2, 0.5],[0.5, 1]]*1e-10 the quadratic p = x^2 - 6xy + y^2 + 10 solves
# div(K grad p) = 0 only because of Kxy, so pp(0.5,0.5) = 9 is recovered exactly
# only if the off-diagonal term is honoured. sat_n is inert (0 everywhere).

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD9
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 1000.0
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

# variable_order = SECOND so the quadratic manufactured solution is in the FE
# space (exact). The Physics applies it to both primary variables; sat_n is held
# Dirichlet 0 on every boundary with a zero IC, so it stays identically 0.
[Physics/VEFlow/FiniteElement]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 2.0e-10
    K_up_yy = 1.0e-10
    K_up_xy = 0.5e-10
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    variable_order = SECOND
    define_geometry_variables = false # z_top/z_bottom declared in [AuxVariables] below
  []
[]

[Functions]
  [p_manufactured]
    type = ParsedFunction
    expression = 'x*x - 6*x*y + y*y + 10'
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 10.0
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
  [pp_all]
    type = FunctionDirichletBC
    variable = pp_top
    boundary = 'left right top bottom'
    function = p_manufactured
  []
  [sat_n_all]
    type = DirichletBC
    variable = sat_n
    boundary = 'left right top bottom'
    value = 0.0
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
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  [pp_center]
    type = PointValue
    variable = pp_top
    point = '0.5 0.5 0'
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
