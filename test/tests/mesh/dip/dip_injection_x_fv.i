# FV version of dip_injection_x.i -- x-dip only, centre-sideset injection.
#
# Finite-volume counterpart: MooseVariableFVReal unknowns, VEFVMassTimeDerivative
# + VEFVAdvectiveFlux, functor relperm. z_top = 0.034899*x (x-dip only); CO2
# injected through the centre 10 m of the left boundary (20 <= y <= 30) migrates
# updip in +x. co2_integral = value*L_inj*t/(H*phi*rho_n) = 1e-5*t (L_inj = 10 m).

[Mesh]
  [base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 50
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 50
  []
  [injection_zone]
    type = ParsedGenerateSideset
    input = base_mesh
    combinatorial_geometry = 'x < 1e-6 & y >= 20.0 & y <= 30.0'
    normal = '-1 0 0'
    new_sideset_name = 'injection_zone'
  []
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

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '0.034899 * x'
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '0.034899 * x - 20'
  []
  [pp_top_init]
    type = ParsedFunction
    expression = '1020 * 9.81 * 0.034899 * (100 - x)'
  []
  [nose_analytical]
    type = ParsedFunction
    expression = '1.0e-12 * 320.0 * 9.81 * 0.034899 * t / (0.2 * 5.0e-5)'
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

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_init
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.0
  []
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = z_top_func
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_func
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
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
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [co2_injection]
    type = FVNeumannBC
    variable = sat_n
    boundary = injection_zone
    value = 2.8e-3
  []
[]

[Materials]
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
  end_time = 5.0e6
  dt = 1.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

[Postprocessors]
  [gravity_nose_x]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_25m]
    type = PointValue
    variable = sat_n
    point = '25 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_at_50m]
    type = PointValue
    variable = sat_n
    point = '50 25 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mass conservation: value * L_inj * t / (H * phi * rho_n), L_inj = 10 m.
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_mass]
    type = FunctionValuePostprocessor
    function = '2.8e-3 * 10.0 * t / (20.0 * 0.2 * 700.0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
