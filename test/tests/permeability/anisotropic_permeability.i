# Anisotropic / off-diagonal permeability verification (method of manufactured solution).
#
# VEPermeability assembles a full symmetric 2x2 tensor ve_K_up = [[Kxx, Kxy],[Kxy, Kyy]],
# but every other test uses isotropic Kxx=Kyy, Kxy=0, so the off-diagonal coupling in the
# Darcy flux (K_up . grad Phi) is otherwise untested. Real Petrel fields carry rotated
# principal axes (Kxy != 0); CLAUDE.md subtlety 6 flags this as must-honour.
#
# Single-phase setup: sat_n = 0 everywhere (kr_n = 0, no CO2 flux; the CO2 equation is
# inert), so the brine pressure equation reduces to the constant-coefficient elliptic
# operator  div(K_up . grad pp) = 0  (flat formation -> no buoyancy; H, kr_w, rho/mu all
# constant factor out).
#
# Manufactured solution: with K_up = [[2, 0.5],[0.5, 1]] * 1e-10,
#   p(x,y) = x^2 - 6 x y + y^2 + 10
# satisfies div(K grad p) = Kxx*p_xx + 2*Kxy*p_xy + Kyy*p_yy
#                         = 2*2 + 2*(0.5)*(-6) + 1*2 = 0   (in units of 1e-10).
# Imposing p on all boundaries, the interior solution is exactly p (second-order LAGRANGE
# reproduces a quadratic exactly under Galerkin, on any mesh). The interior node value
#   p(0.5, 0.5) = 0.25 - 1.5 + 0.25 + 10 = 9.0
# is recovered ONLY if the off-diagonal Kxy enters the operator. If Kxy were dropped the
# quadratic is no longer a solution (operator = 2*2 + 1*2 = 6 != 0) and pp(0.5,0.5) /= 9.

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

[GlobalParams]
  VEDictator = ve_dictator
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
  [ve_dictator]
    type = VEDictator
    porous_flow_vars = 'pp_top sat_n'
    ve_flavour = sharp_interface
  []
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Variables]
  # Second order so the quadratic manufactured solution is in the FE space (exact).
  [pp_top]
    order = SECOND
    family = LAGRANGE
  []
  [sat_n]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  [p_manufactured]
    type = ParsedFunction
    expression = 'x*x - 6*x*y + y*y + 10'
  []
[]

[ICs]
  # Start away from the exact solution so the elliptic solve must actually recover it.
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
    K_up_xx = 2.0e-10
    K_up_yy = 1.0e-10
    K_up_xy = 0.5e-10
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
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  # Interior node value -- recovered only if Kxy is used in the flux operator.
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
