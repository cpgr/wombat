# McWhorter-Sunada capillary verification -- 1D counter-current imbibition, FV
#
# FV counterpart of mcwhorter_fe.i. Exercises the VEFVCapPressure functor and the
# VEFVAdvectiveFlux capillary=true path (grad(Pc^up).n via FV reconstruction).
# Same physics, same reference (mcwhorter_reference.csv); see mcwhorter_fe.i for
# the derivation. Parameters MUST stay in sync with mcwhorter_reference.py.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 100
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
    type = MooseVariableFVReal
  []
  [sat_n]
    type = MooseVariableFVReal
  []
[]

[Functions]
  # axis = x makes this a function of the x-coordinate (not time, the default).
  [mcwhorter_ref]
    type = PiecewiseLinear
    data_file = mcwhorter_reference.csv
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
    value = 0.05
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [z_bottom]
    type = MooseVariableFVReal
    initial_condition = -10.0
  []
[]

[FVKernels]
  [co2_storage]
    type = VEFVMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
    z_top = z_top
    z_bottom = z_bottom
  []
  [co2_flux]
    type = VEFVAdvectiveFlux
    variable = sat_n
    fluid_phase = 0
    pp_top = pp_top
    z_top = z_top
    z_bottom = z_bottom
    capillary = true
  []
  [brine_storage]
    type = VEFVMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
    z_top = z_top
    z_bottom = z_bottom
  []
  [brine_flux]
    type = VEFVAdvectiveFlux
    variable = pp_top
    fluid_phase = 1
    pp_top = pp_top
    z_top = z_top
    z_bottom = z_bottom
  []
[]

# Counter-current imbibition: fixed inlet (left), closed far end (right).
[FVBCs]
  [sat_n_inlet]
    type = FVDirichletBC
    variable = sat_n
    boundary = left
    value = 0.7
  []
  [pp_top_left]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
  []
[]

[Materials]
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
[]

[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n = sat_n
  []
  [cap_pressure]
    type = VEFVCapPressure
    sat_n = sat_n
    z_top = z_top
    z_bottom = z_bottom
    S_wr = 0.0
    delta_rho = 300.0
    gravity = '0 0 -9.81'
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
  end_time = 1.0e5
  dt = 500.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
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
    type = ElementValueSampler
    variable = sat_n
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
