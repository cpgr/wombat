# AD-vs-FD Jacobian consistency of the hysteretic relperm path (UO -> VERelPerm
# adapter -> VEAdvectiveFluxS), evaluated on the imbibition branch.
#
# Uniform interior state: sat_n = 0.5, sat_n_max = 0.8 (a frozen aux). Since
# sat_n < sat_n_max the hysteretic CO2 relperm follows the imbibition scanning
# curve, and 0.5 is strictly between the Land residual s_gr = 0.8/(1+3*0.8) = 0.235
# and the turning point 0.8, so dRelativePermeability is smooth there (away from the
# clamp kink at 0/1, the turning-point kink at sat_n_max, and the trapped kink at
# s_gr -- cf. Gotcha A). This checks the scanning-curve derivative seeded by the
# 3-argument relativePermeabilityAD.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 20
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
    type = VERelPermHysteresisUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
    S_gr_max = 0.25
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
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
    value = 0.5
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
  # Frozen turning point > sat_n, forcing the imbibition branch.
  [sat_n_max]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.8
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
    value = 1.0e5
  []
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
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
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
    sat_n_max = sat_n_max
  []
[]

[Executioner]
  type = Transient
  end_time = 1000.0
  dt = 1000.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_max_its = 20
[]
