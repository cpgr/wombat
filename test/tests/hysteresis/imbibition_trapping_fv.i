# Drainage -> imbibition FV demonstration: residual CO2 trapping by hysteresis.
#
# A 1-D VE column is first drained (CO2 injection from left, Phase 1) then
# imbibed (brine injection from left, Phase 2). The hysteretic relperm
# (VERelPermHysteresisUO) locks residual CO2 in cells swept by the brine front,
# confirming that the scanning-curve logic immobilises CO2 in a real transient solve.
#
# Parameters (equal rho and mu eliminate buoyancy -- pure Buckley-Leverett):
#   K = 1e-10 m2, phi = 0.2, H = 10 m, pp_drop = 1e5 Pa over 100 m
#   rho_n = rho_w = 1000 kg/m3, mu_n = mu_w = 1e-3 Pa.s
#   S_wr = 0.1, S_gr_max = 0.2  =>  C = (1-S_wr)/S_gr_max - 1 = 3.5
#   At full drainage (sat_n_max = 1-S_wr = 0.9, eff s_max = 1):
#     s_gr = s_max/(1+C*s_max) = 2/9,  sat_n_trap(full) = (1-S_wr)*s_gr = 0.2
#
# Phase 1 (0 <= t <= T1 = 20000 s):  sat_n = 1-S_wr = 0.9 at left BC (CO2 inlet).
#   Darcy velocity u = K*(dP/dx)/mu = 1e-10*(1e3)/1e-3 = 1e-4 m/s
#   Piston front (linear relperm, M=1) at x = u*T1/phi = 10 m (10 cells).
#
# Phase 2 (T1 < t <= T2 = 40000 s):  sat_n = 0 at left BC (brine inlet).
#   Brine sweeps the drainage zone. The imbibition scanning curve reduces
#   CO2 mobility; cells swept by the brine front retain sat_n_trap = 0.2 of
#   residually trapped CO2.
#
# Verification: trapped_volume is non-zero and grows during Phase 2 while
# co2_volume decreases, confirming hysteresis does trap CO2 in a real flow solve.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 100
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
  [relperm_uo]
    type = VERelPermHysteresisUO
    S_wr = 0.1
    krn_max = 1.0
    krw_max = 1.0
    S_gr_max = 0.2
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
  # Phase 1: CO2-saturated inlet. Phase 2: brine inlet.
  [sat_n_left_func]
    type = ParsedFunction
    expression = 'if(t <= 20000.0, 0.9, 0.0)'
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
  # sat_n_max: IC must equal initial sat_n (here 0) per VESaturationMaxAux docs.
  [sat_n_max]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [sat_n_trap]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
[]

[AuxKernels]
  # sat_n_max must be advanced first so sat_n_trap sees the updated turning point.
  [sat_n_max_aux]
    type = VESaturationMaxAux
    variable = sat_n_max
    sat_n = sat_n
    execute_on = TIMESTEP_END
  []
  [sat_n_trap_aux]
    type = VETrappedSaturationAux
    variable = sat_n_trap
    relperm_uo = relperm_uo
    sat_n = sat_n
    sat_n_max = sat_n_max
    execute_on = TIMESTEP_END
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
    type = FVFunctionDirichletBC
    variable = sat_n
    boundary = left
    function = sat_n_left_func
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
  end_time = 40000.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  [TimeSteppers]
    [adapt]
      type = IterationAdaptiveDT
      dt = 2000.0
      optimal_iterations = 6
      growth_factor = 1.5
      cutback_factor = 0.5
    []
  []

  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  nl_max_its = 25
  l_tol = 1.0e-8
[]

[Postprocessors]
  [co2_volume]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [trapped_volume]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n_trap
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [satn_max_pp]
    type = ElementExtremeValue
    variable = sat_n
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [smax_max_pp]
    type = ElementExtremeValue
    variable = sat_n_max
    value_type = max
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
