# Unit test for VEGeometry gradient calculation on a planar dipping surface.
#
# Top surface is a plane dipping in both x and y:
#   z_top    = -1000 - 0.01*x - 0.005*y   [m]
#   z_bottom = -1100 - 0.01*x - 0.005*y   [m]
#   H        = z_top - z_bottom = 100 m    (uniform everywhere)
#
# Expected values from VEGeometry:
#   ve_grad_z_top_x = dz_top/dx = -0.01  (1% dip in +x)
#   ve_grad_z_top_y = dz_top/dy = -0.005 (0.5% dip in +y)
#   ve_H            = 100 m
#
# Because z_top is a piecewise-linear LAGRANGE FIRST ORDER field, the
# gradient is exact on any mesh for a planar surface.  A 10x10 mesh is
# sufficient.
#
# The solve itself is trivial: sat_n = 0 everywhere, pp_top fixed at 1e7 Pa
# on all boundaries, so no flux and the Newton solver converges in one step.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
[]

[GlobalParams]
  VEDictator = ve_dictator
[]

[UserObjects]
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = sharp_interface
  []
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.2
    krn_max = 0.8
    krw_max = 1.0
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '-1000 - 0.01*x - 0.005*y'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '-1100 - 0.01*x - 0.005*y'
  []
[]

[Variables]
  [pp_top]
    initial_condition = 1.0e7
  []
  [sat_n]
    initial_condition = 0.0
  []
[]

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
  []
  [grad_z_top_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [grad_z_top_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [H_out]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [z_top_aux]
    type = FunctionAux
    variable = z_top
    function = z_top_func
    execute_on = 'INITIAL'
  []
  [z_bottom_aux]
    type = FunctionAux
    variable = z_bottom
    function = z_bottom_func
    execute_on = 'INITIAL'
  []
  [grad_z_top_x_aux]
    type = MaterialRealVectorValueAux
    variable = grad_z_top_x
    property = ve_grad_z_top
    component = 0
    execute_on = 'TIMESTEP_END'
  []
  [grad_z_top_y_aux]
    type = MaterialRealVectorValueAux
    variable = grad_z_top_y
    property = ve_grad_z_top
    component = 1
    execute_on = 'TIMESTEP_END'
  []
  [H_out_aux]
    type = MaterialRealAux
    variable = H_out
    property = ve_H
    execute_on = 'TIMESTEP_END'
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
  [pp_top_bc]
    type = DirichletBC
    variable = pp_top
    boundary = top
    value = 1.0e7
  []
  [pp_bottom_bc]
    type = DirichletBC
    variable = pp_top
    boundary = bottom
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
  [permeability]
    type = VEPermeability
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
  []
  [fluid_props]
    type = VEFluidPropertiesConst
    rho_co2 = 700.0
    rho_brine = 1020.0
    mu_co2 = 5.0e-5
    mu_brine = 5.0e-4
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [plume_reconstruction]
    type = VEPlumeReconstruction
    sat_n = sat_n
    S_wr = 0.2
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
  [grad_z_top_x_avg]
    type = ElementAverageValue
    variable = grad_z_top_x
    execute_on = 'TIMESTEP_END'
  []
  [grad_z_top_y_avg]
    type = ElementAverageValue
    variable = grad_z_top_y
    execute_on = 'TIMESTEP_END'
  []
  [H_avg]
    type = ElementAverageValue
    variable = H_out
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
