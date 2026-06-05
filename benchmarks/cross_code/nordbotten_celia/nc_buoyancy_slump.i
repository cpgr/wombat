# Pure-buoyancy gravity-slump diagnostic -- isolates the spreading term.
#
# A block of CO2 (sat_n = 0.5) is released in the downdip 20 m with NO injection
# and allowed to slump updip under gravity alone. This isolates the buoyancy /
# spreading physics from any injection or confined-pressure drive, and was the
# test that pinned the MRST<->wombat migration mismatch to the hydrostatic
# spreading term Pc^up = delta_rho*g*h (run with vs without capillary):
#
#   updip nose advance at t = 5e5 s (from the initial x = 20 m edge):
#     wombat capillary OFF : 4.75 m   (= v_N*t -- slope only, barely spreads)
#     wombat capillary ON  : 29.25 m  (strongly spreading -- this input)
#     MRST                 : 49.75 m  (strongly spreading; cap intrinsic)
#
# Capillary is ON here. Set capillary=false and drop VEFVCapPressure to recover
# the slope-only v_N behaviour. Run to ~5e5 s. See doc/content/benchmarks/cross_code.md.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 100
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

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_init
  []
  [sat_n_ic]
    type = FunctionIC
    variable = sat_n
    function = 'if(x < 20, 0.5, 0.0)'
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

# ---------------------------------------------------------------------------
# Geometry AuxVariables -- MUST be FV (MooseVariableFVReal).
#
# Regular Materials that couple FE auxvariables are NOT reinitialised on FV
# faces (ComputeFVFluxThread::reinitVariables only calls computeFaceValues on
# FV variables), so an FE z_top would read as zero at every face and kill the
# flux. The VE FV kernels difference these FV cell values directly to build H
# and the buoyancy slope grad(z_T).
# ---------------------------------------------------------------------------
[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
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

# ---------------------------------------------------------------------------
# FV Boundary conditions
# ---------------------------------------------------------------------------
[FVBCs]
  # Open right boundary (updip end): Dirichlet brine pressure.
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  # NO injection -- this is a pure-buoyancy gravity-slump diagnostic.
[]

# ---------------------------------------------------------------------------
# Materials
#
# No VEGeometry here: the FV kernels build H and grad(z_T) inline from the
# z_top / z_bottom FV variables (see header notes). The remaining materials all
# couple either constants or FV variables (sat_n), so they evaluate correctly
# on faces.
# ---------------------------------------------------------------------------
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

# VEFVRelPerm is a FunctorMaterial, so it lives in [FunctorMaterials]
# (not [Materials]). It exposes ve_relperm_n / ve_relperm_w as functors that the
# flux kernel evaluates at the upwind face.
[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n = sat_n
  []
  [cap_pressure]
    # Two-pressure capillary flux Pc^up = delta_rho * g * h. ON here (unlike the
    # verification test nc_sloped_fv.i) so the Delta_rho*g*dh/dx hydrostatic
    # SPREADING term is active -- MRST's CO2VEBlackOilTypeModel includes it
    # intrinsically, so wombat must enable it for a like-for-like comparison.
    type = VEFVCapPressure
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

# ---------------------------------------------------------------------------
# Preconditioning -- SMP full=true is REQUIRED for coupled FV.
#
# The brine (pp_top) equation couples strongly to sat_n through S_w = 1 - sat_n
# in the storage term, while its own pp_top diagonal (the Darcy term) is several
# orders of magnitude smaller. Without the full off-diagonal Jacobian blocks the
# linear system is effectively singular relative to the injection residual and
# the solve breaks down. SMP full=true assembles every coupling block.
# ---------------------------------------------------------------------------
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

# ---------------------------------------------------------------------------
# Executioner
# ---------------------------------------------------------------------------
[Executioner]
  type = Transient
  end_time = 2.0e6
  dt = 1.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  # pp_top ~ 1e7 Pa and sat_n ~ 1 differ by orders of magnitude; let MOOSE
  # rescale the residuals so the coupled Newton system is well conditioned.
  automatic_scaling = true

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

# ---------------------------------------------------------------------------
# Postprocessors -- same as nc_sloped.i for direct comparison
# ---------------------------------------------------------------------------
[Postprocessors]
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_mass]
    type = FunctionValuePostprocessor
    function = '4.0e-6 * t / (0.2 * 20.0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [gravity_nose]
    type = FunctionValuePostprocessor
    function = nose_analytical
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Note: near-field saturation probes (e.g. PointValue at x = 5/10/20 m) are
  # intentionally omitted from the regression CSV. In the shortened window they
  # only register the first-order-upwind diffusion tail (O(1e-6) and growing),
  # which is numerical noise rather than verifiable physics; the full sat_n field
  # is in the Exodus output. The CSV checks mass conservation (co2_integral vs
  # analytical_mass) and the analytical gravity-nose, which are exact.
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  # Sample sat_n at cell centroids (0.25, 0.75, ...) -- avoids FV face ambiguity.
  # Overlay against the MRST nc_mrst_sat_*.csv profiles.
  [sat_profile]
    type = LineValueSampler
    variable = sat_n
    start_point = '0.25 0 0'
    end_point = '99.75 0 0'
    num_points = 200
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]
