# McWhorter-Sunada capillary verification, CAPILLARY-FRINGE flux -- 1D, FE
#
# Fringe counterpart of mcwhorter_fe.i. Same 1D counter-current capillary
# imbibition (constant densities, flat formation, closed far boundary => q_t = 0)
# so the CO2 saturation obeys phi dS/dt = d/dx[D(S) dS/dx]. The capillary slope
# that enters D(S) is the FRINGE one:
#
#   D(S) = (K/mu) S (1-S) * (rho_w-rho_n) g H / S_n(h(S)),
#   S_n(h) = 1 - Sw(delta_rho g h),  h(S) the table column inversion.
#
# This is D_sharp(S)/S_n(h(S)) -- the column imbibes FASTER than sharp, a clearly
# different profile, so it verifies the fringe S_n(h) factor's SIGN and MAGNITUDE
# (which neither the Jacobian tests nor the linear-table no-flow equilibrium can).
# Reference: the Boltzmann self-similar profile from mcwhorter_fringe_reference.py
# (cross-checked against an independent method-of-lines solve to L2 ~ 1e-7).
#
# H = 20 m keeps S_n(h(S)) ~ 0.6..0.7 over the active band (away from the
# D ~ sqrt(S) degeneracy at S -> 0). co2_l2_error is ||sat_n - S_ref|| at
# t = 4e4 s; gold 0 with abs_zero = 0.012 asserts a match (the SHARP coefficient
# would give L2 ~ 0.022 -> fail; a wrong SIGN gives anti-diffusion -> blow-up).
#
# Parameters MUST stay in sync with mcwhorter_fringe_reference.py.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 400
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
  [ve_pc_table]
    # MUST match PC_POINTS / SW_POINTS in mcwhorter_fringe_reference.py.
    type = VECapillaryPressureTableUO
    sat_lr = 0.2
    pc_points = '0  15696  31392  47088  62784'
    sw_points = '1.0  0.8  0.6  0.4  0.2'
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
  []
[]

[Functions]
  [mcwhorter_ref]
    type = PiecewiseLinear
    data_file = mcwhorter_fringe_reference.csv
    format = columns
    axis = x
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
    value = 0.04 # Si
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
    initial_condition = -20.0 # H = 20 m
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
    capillary = true
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
  [sat_n_inlet]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 0.3 # S0
  []
  [pp_top_left]
    type = DirichletBC
    variable = pp_top
    boundary = left
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
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
  []
  [saturation]
    type = VESaturation
    sat_n = sat_n
  []
  [plume_reconstruction]
    type = VEPlumeReconstruction
    mode = capillary_fringe
    sat_n = sat_n
    S_wr = 0.0
    pc_uo = ve_pc_table
  []
  [cap_pressure]
    type = VEUpscaledCapPressure
    mode = capillary_fringe
    sat_n = sat_n
    pc_uo = ve_pc_table
    gravity = '0 0 -9.81'
    S_wr = 0.0
  []
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
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
  end_time = 4.0e4
  dt = 200.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-9
  nl_abs_tol = 1.0e-11
  nl_max_its = 20
  l_tol = 1.0e-8
[]

[Postprocessors]
  [co2_l2_error]
    type = ElementL2Error
    variable = sat_n
    function = mcwhorter_ref
    execute_on = 'FINAL'
  []
[]

[VectorPostprocessors]
  [saturation_profile]
    type = LineValueSampler
    variable = sat_n
    start_point = '0 0 0'
    end_point = '100 0 0'
    num_points = 401
    sort_by = x
    execute_on = 'FINAL'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
