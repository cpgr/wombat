# Cross-code benchmark, Tier A1 -- conceptual anticline, stripped physics. FV PATH.
#
# Reads the shared 2D top-surface grid written by anticline_grid_generator.i and runs
# the VE flow solve whose observables are compared against MRST co2lab (see
# doc/content/benchmarks/cross_code.md). Tier A1 is the workhorse: it exercises the
# sloped-surface buoyancy drive in genuine 2D plus a structural trap, with constant
# densities (no EOS ambiguity) and NO trapping/dissolution.
#
# This is the FINITE-VOLUME version (VEFVMassTimeDerivative + VEFVAdvectiveFlux,
# functor relperm VEFVRelPerm). The FV path is the robust choice for the
# advection-dominated migration the full benchmark needs -- FE continuous Galerkin
# oscillates here (HANDOVER, hysteresis Phase 4b note 2). Mirrors test/tests/mesh/dip/
# dip_injection_fv.i, but with the geometry/petrophysics read from the shared Exodus
# grid instead of analytic functions.
#
# Physics (must match co2lab exactly):
#   sharp interface, S_wr = 0.2
#   rho_CO2 = 700, rho_brine = 1000 kg/m^3   mu_CO2 = 6e-5, mu_brine = 8e-4 Pa.s
#   capillary ON (two-pressure spreading term, to match MRST), dissolution off, trapping off
#
# Geometry: +x is updip; CO2 injected at the downdip (x = 0) edge migrates +x and
# pools under the dome crest near x = 7 km.
#
# FIELD PLUMBING -- two FV gotchas govern the choices below (full writeup, including
# why each one makes the matrix singular, in doc/content/benchmarks/cross_code.md,
# "FV input gotchas"):
#   1. Fields a flux kernel touches on a FACE must be MooseVariableFVReal, not FE
#      auxvars (which read zero on faces). So phi_bar/K_up are read from Exodus into FV
#      variables. VEGeometry is the inverse case: it feeds only the elemental mass PP
#      (ve_H) and MUST couple FE geometry, because its coupledGradient is empty for an
#      FV var (out-of-bounds crash). The face rule does not apply to an elemental
#      consumer, so FE geometry there is correct.
#   2. FV functor geometry must be live before the first Jacobian -> set z_top_fv/
#      z_bottom_fv by FunctionIC (analytic Tier A surface, identical to the generator),
#      not an INITIAL AuxKernel. Tier B/C will use initial_from_file_var (also an IC).
#
# STATUS: scaffold. Schedule is now physical: 1 Mt/yr through a narrow 400 m inlet for
# 10 yr (shut off by [Controls]), then ~190 yr of free updip migration into the dome.
# PRODUCTION still needs:
#   - extend end_time toward 500 yr to let the plume fully pool/trap in the dome,
#   - consider a TVD limiter (advected_interp_method = vanLeer) for the sharp front,
#   - compare (CSVDiff) against reference/anticline_reference_*.csv from co2lab.
# Tune solver/timestepping live (the user runs ./run_tests).

[Mesh]
  [base]
    type = FileMeshGenerator
    file = anticline_grid_generator_out.e
    use_for_exodus_restart = true
  []
  # CO2 injection patch: 400 m of the downdip (x = 0) edge, centred at y = 5 km.
  # Narrow inlet -> the fixed 1 Mt/yr enters at high local saturation and migrates
  # updip as a long, focused finger rather than a broad sheet.
  [injection_zone]
    type = ParsedGenerateSideset
    input = base
    combinatorial_geometry = 'x < 1e-3 & y >= 4800 & y <= 5200'
    normal = '-1 0 0'
    new_sideset_name = 'injection_zone'
  []
[]
[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 6.0e-5
  []
  [brine_fp]
    type = ConstantFluidProperties
    density = 1000.0
    viscosity = 8.0e-4
  []
[]

[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    # S_wr = 0 here gives the SIMPLE (non-renormalised) sharp curves kr_n = sat_n,
    # kr_w = 1 - sat_n, matching MRST co2lab's 'sharp_interface_simple' relperm
    # (krG = sG, krW = sw). The renormalised form (S_wr = 0.2 -> kr_n = sat_n/0.8)
    # is more accurate but drives kr_w -> 0 at a filled column, which stalls MRST's
    # stock solver; the simple form keeps kr_w >= 0.2 and lets both codes complete.
    # The physical residual (0.2) is still applied to the column height via the
    # capillary material [FunctorMaterials/cap_pressure] (S_wr = 0.2).
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

[AuxVariables]
  # --- Petrophysics read from the shared Exodus grid as FV variables (so they
  #     evaluate correctly on FV faces -- see gotcha 1 in the header). The generator
  #     writes these elemental, a clean elemental -> FV restart. ---
  [phi_bar]
    type = MooseVariableFVReal
    initial_from_file_var = phi_bar
    initial_from_file_timestep = LATEST
  []
  [K_up_xx]
    type = MooseVariableFVReal
    initial_from_file_var = K_up_xx
    initial_from_file_timestep = LATEST
  []
  [K_up_yy]
    type = MooseVariableFVReal
    initial_from_file_var = K_up_yy
    initial_from_file_timestep = LATEST
  []
  [K_up_xy]
    type = MooseVariableFVReal
    initial_from_file_var = K_up_xy
    initial_from_file_timestep = LATEST
  []

  # --- FV functor geometry for the FV flux kernels (set by FunctionIC below). ---
  [z_top_fv]
    type = MooseVariableFVReal
  []
  [z_bottom_fv]
    type = MooseVariableFVReal
  []

  # --- FE (LAGRANGE) geometry read from Exodus, ONLY for VEGeometry -> ve_H (the
  #     elemental mass PP). VEGeometry calls coupledGradient, which is empty for an FV
  #     variable (out-of-bounds crash), so it needs gradient-capable FE geometry. Safe
  #     because ve_H is never evaluated on a face -- see gotcha 1 in the header.
  #     VEGeometry now needs both surfaces (ve_H = z_top - z_bottom); the mesh carries
  #     nodal z_top and H, so z_bottom is derived as z_top - H by the AuxKernel below. ---
  [z_top]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = z_top
    initial_from_file_timestep = LATEST
  []
  [H]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = H
    initial_from_file_timestep = LATEST
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Nodal z_bottom = z_top - H for VEGeometry. Static geometry, so INITIAL is enough;
  # reproduces the exact ve_H = H and ve_grad_H = grad(z_top) - grad(z_bottom) = grad(H)
  # that VEGeometry's removed thickness mode read from the H field directly.
  [z_bottom_aux]
    type = ParsedAux
    variable = z_bottom
    coupled_variables = 'z_top H'
    expression = 'z_top - H'
    execute_on = 'INITIAL'
  []
[]

[Functions]
  # Analytic anticline surface -- IDENTICAL to anticline_grid_generator.i. Used to set
  # the FV functor geometry as an initial condition (live before the first assembly).
  [ztop_fn]
    type = ParsedFunction
    expression = '-800 + 0.0261859*x + 40*exp(-((x-7000)^2 + (y-5000)^2)/(2*1200^2))'
  []
  [zbottom_fn]
    type = ParsedFunction
    expression = '-850 + 0.0261859*x + 40*exp(-((x-7000)^2 + (y-5000)^2)/(2*1200^2))'
  []
  # Hydrostatic brine pressure referenced to the open updip boundary. Approximate
  # (ignores the dome) -- it is only an IC; the brine solve relaxes it.
  [pp_top_init]
    type = ParsedFunction
    expression = '1000 * 9.81 * 0.0261859 * (10000 - x)'
  []
[]

[ICs]
  # FV geometry as initial conditions (closed/ghosted before assembly -- see header).
  [z_top_fv_ic]
    type = FunctionIC
    variable = z_top_fv
    function = ztop_fn
  []
  [z_bottom_fv_ic]
    type = FunctionIC
    variable = z_bottom_fv
    function = zbottom_fn
  []
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
[]

[FVKernels]
  # CO2 mass equation (variable: sat_n, phase 0).
  [co2_storage]
    type = VEFVMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
    z_top = z_top_fv
    z_bottom = z_bottom_fv
  []
  [co2_flux]
    type = VEFVAdvectiveFlux
    variable = sat_n
    fluid_phase = 0
    pp_top = pp_top
    z_top = z_top_fv
    z_bottom = z_bottom_fv
    # Two-pressure capillary flux ON: activates the hydrostatic spreading term
    # delta_rho*g*grad(h). MRST's CO2VEBlackOilTypeModel includes this
    # intrinsically, so it must be on here for a like-for-like comparison --
    # it controls how the plume spreads, fills, and spills the dome trap.
    # (See doc/content/benchmarks/cross_code.md and the nordbotten_celia inputs.)
    capillary = true
  []
  # Brine mass equation (variable: pp_top, phase 1).
  [brine_storage]
    type = VEFVMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
    z_top = z_top_fv
    z_bottom = z_bottom_fv
  []
  [brine_flux]
    type = VEFVAdvectiveFlux
    variable = pp_top
    fluid_phase = 1
    pp_top = pp_top
    z_top = z_top_fv
    z_bottom = z_bottom_fv
  []
[]

[FVBCs]
  # Open updip boundary (right): brine exits via the Dirichlet pressure; CO2 exits via
  # the advective-outflow BC. Without the CO2 outflow, any plume that spills updip past
  # the dome hits a no-flow wall and smears laterally along x = 10 km -- this lets it
  # leave the domain at the Darcy rate instead.
  [pp_updip]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 0.0
  []
  [co2_outflow]
    type = VEFVAdvectiveOutflowBC
    variable = sat_n
    boundary = right
    fluid_phase = 0
    pp_top = pp_top
    z_top = z_top_fv
    z_bottom = z_bottom_fv
  []
  # CO2 injection on the downdip patch -- 1 Mt/yr (typical industrial CCS rate),
  # switched OFF after 10 yr by the [Controls] block so the plume migrates freely
  # afterwards. The CO2 residual is in mass units (VEFVMassTimeDerivative storage is
  # H*phi*rho*S), so this FVNeumannBC value is a mass flux per unit injection-edge
  # length [kg/(m*s)]; value * edge_length = total mass rate. The injection_zone edge
  # is now 400 m (y in [4800,5200]), so to hold the SAME 1 Mt/yr total:
  #   value = (1e9 kg / 3.1557e7 s) / 400 m = 31.69 kg/s / 400 m = 7.922e-2 kg/(m*s).
  # Match this to co2lab's rate well when wiring the cross-code comparison.
  [co2_injection]
    type = FVNeumannBC
    variable = sat_n
    boundary = injection_zone
    value = 7.922e-2
  []
[]

# Injection window: enable co2_injection only for [0, 10 yr], then disable it so the
# downdip boundary reverts to natural zero-flux and the plume migrates updip into the
# dome. set_sync_times lands a timestep exactly on the 10 yr shut-off, so the total
# injected mass is clean. (reverse_on_false defaults true -> disabled after the window.)
[Controls]
  [injection_window]
    type = TimePeriod
    enable_objects = 'FVBCs/co2_injection'
    implicit = false
    start_time = 0.0
    end_time = 3.1557e8 # 10 yr
    set_sync_times = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Materials]
  # ve_H for the (elemental) mobile-mass PP. Couples the FE geometry: VEGeometry takes
  # coupledGradient(z_top), empty for an FV var. Consumed only at element centroids, so
  # the FV face-value rule does not apply (see header gotcha 1). z_bottom is the derived
  # nodal field z_top - H (see [AuxKernels/z_bottom_aux]).
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
  []
  [porosity]
    type = VEPorosity
    phi_bar = phi_bar
  []
  [permeability]
    type = VEPermeability
    K_up_xx = K_up_xx
    K_up_yy = K_up_yy
    K_up_xy = K_up_xy
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
  # Provides ve_dPcup_dsatn / ve_dPcup_dH for the capillary=true CO2 flux above.
  # FV functor geometry (z_top_fv/z_bottom_fv) since it is evaluated on faces;
  # S_wr = 0.2 matches the VERelPermSharpUO so the sharp-interface column height
  # h = sat_n*H/(1-S_wr) is consistent between the relperm and the spreading term.
  [cap_pressure]
    type = VEFVCapPressure
    sat_n = sat_n
    z_top = z_top_fv
    z_bottom = z_bottom_fv
    S_wr = 0.2
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
  # 10 yr injection + ~190 yr post-injection migration (~200 yr total). At the
  # gravity-nose speed v_N = K*drho*g*sin(theta)/(phi*mu_n) ~ 1.3e-6 m/s this lets the
  # plume reach the dome crest near x = 7 km. Extend toward 500 yr to watch it pool.
  end_time = 6.3e9

  # Adaptive stepping: small steps while injecting / at the shut-off transition, then
  # grow through the slow migration phase. dtmax keeps the migrating front resolved.
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e6
    growth_factor = 1.4
    cutback_factor = 0.5
    optimal_iterations = 8
  []
  dtmax = 1.0e8

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

# ---------------------------------------------------------------------------
# Benchmark observables (compared against co2lab)
# ---------------------------------------------------------------------------
[Postprocessors]
  # Total mobile CO2 mass [kg] -- primary cross-code metric.
  [mobile_co2_mass]
    type = VEMobileCO2MassPostprocessor
    sat_n = sat_n
    rho_co2 = 700.0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Plume footprint area [m^2].
  [footprint]
    type = VEPlumeFootprintPostprocessor
    sat_n = sat_n
    threshold = 1.0e-3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Mass-conservation cross-check (domain integral of sat_n).
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# Updip nose / thickness profile along the injection centreline.
# sat_n is a discontinuous (FV) variable, so sampling on an element FACE is ambiguous
# (the LineValueSampler warning). On the 50x50 / 10 km mesh (dx = dy = 200), element
# centroids sit at 100, 300, ..., 9900. So sample a centroid row (y = 4900, within the
# injection band [4800,5200]) with 50 points at x = 100..9900 -- every point is an
# element interior, no faces touched.
[VectorPostprocessors]
  [centreline_sat]
    type = LineValueSampler
    variable = sat_n
    start_point = '100 4900 0'
    end_point = '9900 4900 0'
    num_points = 50
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
