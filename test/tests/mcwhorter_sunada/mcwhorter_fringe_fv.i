# McWhorter-Sunada capillary verification, CAPILLARY-FRINGE flux -- 1D, FV
#
# FV counterpart of mcwhorter_fringe_fe.i: exercises the VEFVCapPressure
# capillary_fringe path (per-face Newton inversion + S_n(h) coefficients) in the
# VEFVAdvectiveFlux capillary=true flux. Same physics, same self-similar reference
# (mcwhorter_fringe_reference.csv). See mcwhorter_fringe_fe.i for the derivation.
# The phase velocities are non-zero (counter-current imbibition), so the default
# upwind interpolation is differentiable here -- no average needed.
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
    type = VECapillaryPressureTableUO
    sat_lr = 0.2
    pc_points = '0  15696  31392  47088  62784'
    sw_points = '1.0  0.8  0.6  0.4  0.2'
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
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [z_bottom]
    type = MooseVariableFVReal
    initial_condition = -20.0 # H = 20 m
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

[FVBCs]
  [sat_n_inlet]
    type = FVDirichletBC
    variable = sat_n
    boundary = left
    value = 0.3 # S0
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
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
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
    mode = capillary_fringe
    pc_uo = ve_pc_table
    sat_n = sat_n
    z_top = z_top
    z_bottom = z_bottom
    S_wr = 0.0
    gravity = '0 0 -9.81'
  []
  [fv_fluid_props]
    type = VEFVFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
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
  line_search = basic
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
