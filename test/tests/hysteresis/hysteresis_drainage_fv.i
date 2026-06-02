# Hysteretic relperm in a real FV transient solve: monotonic drainage reduces to
# the sharp Buckley-Leverett result. This exercises the FV hysteresis adapter path
# (VERelPermHysteresisUO -> VEFVRelPerm -> VEFVAdvectiveFlux with the sat_n_max
# functor) inside an actual flow solve, complementing the algebraic column_cycle
# (solve=false) and the FE Jacobian (hysteresis_jacobian) checks.
#
# CO2 is injected into a brine-filled column (sat_n: 0 -> 1, pure drainage). Because
# sat_n only ever rises, sat_n >= sat_n_max everywhere and VERelPermHysteresisUO stays
# on its drainage branch, where it is identical to VERelPermSharpUO. The conserved
# integral co2_volume = integral(sat_n) therefore equals the analytical front 5e-4*t,
# exactly as in buckley_leverett/bl_flat_fv.i -- confirming the hysteretic UO does not
# perturb drainage and that the sat_n_max wiring is correct in the FV functor adapter.
#
# (The imbibition/trapping branch is verified algebraically in column_cycle; a full
# imbibition-displacement flow solve is a separate, harder problem -- CO2 banking to
# the maximum saturation drives kr_w -> 0 and locks the displacement -- and is left as
# a future verification once anti-banking measures, e.g. S_wr > 0 and TVD limiting,
# are in place.)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
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
    ve_flavour = with_hysteresis
  []
  [relperm_uo]
    type = VERelPermHysteresisUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
    S_gr_max = 0.3
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
  [pp_top_linear]
    type = ParsedFunction
    expression = '1.0e5 * (1.0 - x / 100.0)'
  []
  [front_analytical]
    type = ParsedFunction
    expression = '5.0e-4 * t'
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
    value = 0.0
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
  [sat_n_max]
    type = MooseVariableFVReal
    # Drainage from sat_n = 0, so the turning point also starts at 0.
    initial_condition = 0.0
  []
[]

[AuxKernels]
  [sat_n_max_aux]
    type = VESaturationMaxAux
    variable = sat_n_max
    sat_n = sat_n
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
    value = 1.0
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
    sat_n_max = sat_n_max
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
  end_time = 5000.0
  dt = 1000.0

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
  [co2_volume]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_front]
    type = FunctionValuePostprocessor
    function = front_analytical
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
