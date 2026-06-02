# Buckley-Leverett verification test -- FV formulation, 1D flat formation, M = 1
#
# FV counterpart of bl_flat.i. Used as a convergence-bisection diagnostic for
# the VE FV kernels:
#   * Flat formation (z_top = 0, z_bottom = -10) => grad(z_T) = 0, so the
#     buoyancy / surface-slope term is identically zero. If this test converges
#     but nc_sloped_fv.i does not, the problem is in the slope/buoyancy term.
#   * Dirichlet sat_n = 1 inlet (not a Neumann injection) and Dirichlet pressure
#     on BOTH ends -- a much better-conditioned problem than nc_sloped_fv.i. If
#     this ALSO fails to converge, the problem is fundamental to the FV storage
#     or flux Jacobian, independent of geometry and injection.
#
# M = 1 (equal density and viscosity, linear relperm) => f_n = S_n, wave speed:
#   v = K_up * |dp/dx| / (phi * mu)
#     = 1e-10 * 1000 / (0.2 * 1e-3) = 5e-4 m/s
# At t = 1e5 s the front is at x = 50 m, so integral(sat_n, 0, 100) = 50.

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
    ve_flavour = sharp_interface
  []
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
  []
[]

# ---------------------------------------------------------------------------
# Primary variables -- FV cell-average unknowns
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Geometry AuxVariables -- FV, constant (flat formation, no topography).
# grad(z_T) = 0 because z_top is uniform, so the buoyancy term vanishes.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# FV Kernels
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# FV Boundary conditions
# ---------------------------------------------------------------------------
[FVBCs]
  # Pressure: Dirichlet on both ends to fix the pressure gradient.
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
  # Saturation: Dirichlet at the inlet (x = 0). The right boundary is left
  # natural (closed); the front does not reach x = 100 within the test window.
  [sat_n_left]
    type = FVDirichletBC
    variable = sat_n
    boundary = left
    value = 1.0
  []
[]

# ---------------------------------------------------------------------------
# Materials -- no VEGeometry; the FV kernels build H inline from z_top/z_bottom.
# ---------------------------------------------------------------------------
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

# VEFVRelPerm is a FunctorMaterial, so it lives in [FunctorMaterials]
# (not [Materials]). It exposes ve_relperm_n / ve_relperm_w as functors that the
# flux kernel evaluates at the upwind face.
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
  end_time = 1.0e5
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

# ElementValueSampler (not LineValueSampler): sat_n is a cell-constant FV
# variable, discontinuous across faces. Sampling at element centroids avoids the
# "discontinuous variable sampled on a face" warning and gives the cell values
# directly.
[VectorPostprocessors]
  [saturation_profile]
    type = ElementValueSampler
    variable = sat_n
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
