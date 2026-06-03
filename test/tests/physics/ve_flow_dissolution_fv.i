# VE flow [Physics] -- FV, enable_dissolution = true.
#
# Physics-generated counterpart of dissolution/constant_rate_fv.i. enable_dissolution =
# true creates the full dissolution chain (VEDissolution material, c_diss aux variable,
# VEDissolvedCO2Aux accumulator, VEFVDissolutionSink) plus the geometry aux variables;
# only ICs/BCs remain user-written. CO2 is immobile (krn_max = 0) so the CO2 equation
# reduces to the storage + dissolution ODE. CSVDiff against constant_rate_fv.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
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
    krn_max = 0.0 # CO2 immobile: isolates storage + dissolution sink (no advection)
    krw_max = 1.0
  []
[]

[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-10
    K_up_yy = 1.0e-10
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    enable_dissolution = true
    dissolution_flux = 1.4
    dissolution_s_ref = 0.05
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
  [z_top_ic]
    type = FunctionIC
    variable = z_top
    function = '0.0'
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = '-10.0'
  []
[]

[FVBCs]
  [pp_left]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
  []
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 1.0e6
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
  end_time = 100.0
  dt = 20.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-9
  nl_abs_tol = 1.0e-11
  nl_max_its = 20
[]

[Postprocessors]
  [sat_n_avg]
    type = ElementAverageValue
    variable = sat_n
    execute_on = 'FINAL'
  []
  [dissolved_co2_mass]
    type = ElementIntegralVariablePostprocessor
    variable = c_diss
    execute_on = 'FINAL'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
