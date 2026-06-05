# Verification of harmonic-mean K face averaging in VEFVAdvectiveFlux.
#
# Two-cell 1-D mesh with a K jump at x = 1:
#   K1 = 1e-10 m^2 (left cell, x in [0,1])
#   K2 = 1e-12 m^2 (right cell, x in [1,2])
#
# Flat domain (z_top = 0, z_bottom = -10), gravity = 0, sat_n = 0 (brine only).
# With constant density/viscosity and kr_w = 1 (sat_n = 0), the steady-state brine
# equation reduces to div(K grad p) = 0 -- the classical two-point transmissibility
# problem.
#
# TPFA harmonic K at the interior face (equal half-distances d = 0.5 m):
#   K_h = 2*K1*K2/(K1+K2) = 2e-22/1.01e-10 = 1.9802e-12 m^2
#
# Analytic steady-state pressures at the cell centroids (p_L = 1e5, p_R = 0):
#   p1 = 1e5 * 201/202 = 99504.950495... Pa
#   p2 = 1e5 * 100/202 = 49504.950495... Pa
#
# With ConstantFluidProperties the storage term (rho * phi * H * dS/dt) is zero
# because sat_n = 0 is fixed by the Dirichlet BCs, so a single transient step
# from the pp_top = 0 IC converges exactly to the steady-state TPFA solution in
# one Newton iteration (linear problem).
#
# The K_up_xx / K_up_yy AuxVariables MUST be MooseVariableFVReal (not plain CONSTANT
# MONOMIAL) so that VEPermeability's coupledValue() is reinitialised on the neighbour
# element when VEFVAdvectiveFlux calls getNeighborMaterialProperty. A plain FE aux
# variable reads as zero on the neighbour side and collapses the harmonic average to
# zero; the Physics action now errors explicitly if a non-FV variable is supplied.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
  xmin = 0
  xmax = 2
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 5.0e-5
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1020.0
    viscosity = 8.0e-4
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

# K_up_xx / K_up_yy must be MooseVariableFVReal so MOOSE reinitialises them on
# neighbour elements when VEFVAdvectiveFlux calls getNeighborMaterialProperty.
[AuxVariables]
  [K_up_xx]
    type = MooseVariableFVReal
    family = MONOMIAL
    order = CONSTANT
  []
  [K_up_yy]
    type = MooseVariableFVReal
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [K_fn]
    type = ParsedFunction
    expression = 'if(x < 1, 1e-10, 1e-12)'
  []
[]

[AuxKernels]
  [K_up_xx_aux]
    type = FunctionAux
    variable = K_up_xx
    function = K_fn
    execute_on = INITIAL
  []
  [K_up_yy_aux]
    type = FunctionAux
    variable = K_up_yy
    function = K_fn
    execute_on = INITIAL
  []
[]

# Action declares z_top/z_bottom as MooseVariableFVReal; ICs set their values.
[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = K_up_xx
    K_up_yy = K_up_yy
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    gravity = '0 0 0'
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 0.0
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.0
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
    value = 0.0
  []
  [sat_n_right]
    type = FVDirichletBC
    variable = sat_n
    boundary = right
    value = 0.0
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
  dt = 1.0
  num_steps = 1

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1.0e-12
  nl_abs_tol = 1.0e-12
[]

[Postprocessors]
  [pp_cell1]
    type = PointValue
    variable = pp_top
    point = '0.5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pp_cell2]
    type = PointValue
    variable = pp_top
    point = '1.5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = TIMESTEP_END
  []
[]
