# VEGeometry thickness-mode verification (the real-field / Petrel-Exodus path).
#
# Identical to buckley_leverett/bl_flat.i EXCEPT the formation thickness is supplied
# directly as the coupled AuxVariable H (thickness mode) instead of being derived from
# z_top - z_bottom (top_bottom mode). With H = 10 m this is geometrically identical to
# bl_flat (z_top = 0, z_bottom = -10), so the Buckley-Leverett mass balance must be
# byte-identical: co2_volume tracks the same FE-smeared front as bl_flat_csv.csv.
#
# This exercises the thickness-mode branch of VEGeometry (ve_H = H via coupledValue("H"),
# ve_grad_H via coupledGradient("H")), which no other test covers and which real upscaled
# meshes use (H is a native nodal field alongside phi_bar, K_up).

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
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  # Thickness mode: H supplied directly (no z_bottom). LAGRANGE FIRST so grad(H) is
  # available (here uniform -> grad(H) = 0).
  [H]
    order = FIRST
    family = LAGRANGE
    initial_condition = 10.0
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
  [sat_n_left]
    type = DirichletBC
    variable = sat_n
    boundary = left
    value = 1.0
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    H = H
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
  []
[]

[Executioner]
  type = Transient
  end_time = 5000.0
  dt = 1000.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-12
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
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
