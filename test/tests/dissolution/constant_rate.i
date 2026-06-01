# Constant-rate dissolution verification (mass conservation + analytic decay).
#
# Isolates the dissolution machinery (VEDissolution -> VEDissolutionSink -> VEDissolvedCO2Aux)
# from two-phase flow: a single-element mesh with pp_top held Dirichlet at both nodes leaves
# no free pressure DOF (brine equation fully constrained) and grad(pp) = 0, so the advective
# flux vanishes and the CO2 mass equation collapses to a pure ODE:
#
#   d/dt(H phi rho_c sat_n) = -q0 * gate(sat_n)
#
# With H=10, phi=0.2, rho_c=700 (H phi rho_c = 1400), q0 = 1.4 kg/m^2/s, sat_n0 = 0.5, and
# s_ref = 0.05 (so gate = 1 throughout, sat_n stays >> s_ref):
#
#   d(sat_n)/dt = -q0/1400 = -1e-3 /s   ->   sat_n(100 s) = 0.5 - 0.1 = 0.4   (BE-exact:
#                                            constant RHS, so backward Euler is exact)
#   c_diss(t) = q0 * t  = 140 kg/m^2    ->   dissolved mass (area = 1) = 140 kg
#   mobile CO2 mass = rho_c phi H sat_n = 1400 * 0.4 = 560 kg
#   conservation: 560 (mobile) + 140 (dissolved) = 700 = initial free CO2.  (mass conserved)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
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

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e6
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
  # Elemental: VEDissolvedCO2Aux reads the ve_dissolution_rate material property, which is
  # defined at quadrature points (not nodes), so c_diss must be a constant-monomial aux.
  [c_diss]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
[]

[AuxKernels]
  [c_diss_aux]
    type = VEDissolvedCO2Aux
    variable = c_diss
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
  [co2_dissolution]
    type = VEDissolutionSink
    variable = sat_n
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
  # Fix pressure at both nodes -> no free pp DOF (brine equation constrained), grad(pp)=0.
  [pp_left]
    type = DirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
  []
  [pp_right]
    type = DirichletBC
    variable = pp_top
    boundary = right
    value = 1.0e6
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
    type = VEFluidPropertiesConst
    rho_co2 = 700.0
    rho_brine = 1000.0
    mu_co2 = 1.0e-3
    mu_brine = 1.0e-3
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
  [dissolution]
    type = VEDissolution
    dissolution_flux = 1.4
    s_ref = 0.05
  []
[]

[Executioner]
  type = Transient
  end_time = 100.0
  dt = 20.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
  nl_max_its = 20
[]

[Postprocessors]
  [sat_n_avg]
    type = ElementAverageValue
    variable = sat_n
    execute_on = 'timestep_end final'
  []
  [dissolved_co2_mass]
    type = ElementIntegralVariablePostprocessor
    variable = c_diss
    execute_on = 'timestep_end final'
  []
  [mobile_co2_mass]
    type = VEMobileCO2MassPostprocessor
    sat_n = sat_n
    rho_co2 = 700.0
    execute_on = 'timestep_end final'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
