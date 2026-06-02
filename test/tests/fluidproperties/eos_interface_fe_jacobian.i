# Jacobian verification of the FE solve with eos_reference_depth = interface.
#
# In interface mode VEFluidProperties evaluates the EOS at the CO2-brine contact
# z_T + h, with the hydrostatic pressure p = pp_top + rho_n(pp_top)*|g|*h and the
# sharp-interface thickness h = sat_n*H/(1-S_wr). This makes ve_density and
# ve_viscosity depend on BOTH pp_top (through the top-surface pressure and the
# increment density) AND sat_n (through h) -- a coupling that the top_surface mode
# does not have. PetscJacobianTester confirms that both derivatives are seeded
# correctly through the AD chain into the storage and flux kernels.
#
# The state is held at a strictly-interior uniform sat_n = 0.5 (off the sharp-relperm
# kink, cf. HANDOVER Gotcha A) with pressure-dependent fluids and a pressure gradient,
# so every density/viscosity and the interface increment are live and differentiable.
#
# The fluids are SimpleFluidProperties (rho = density0*exp(p/bulk_modulus - ...)), which
# has EXACT analytic derivatives, NOT CO2FluidProperties. This is deliberate. In interface
# mode the density derivative is a composition with a nested coefficient:
#   d(rho)/d(pp_top) = drho/dp|_p * (1 + |g|*h * drho/dp|_{p_top}),
# so any inconsistency d between an EOS's value and its analytic drho/dp is AMPLIFIED by
# |g|*h (~100 here). CO2FluidProperties inverts Span-Wagner with an internal Newton solve;
# its d is tiny (passes top_surface Jacobian tests, where the derivative is just drho/dp),
# but |g|*h * d lands ~1e-3 -- and once that density sits in the heavily-weighted
# mass-STORAGE term it pushes the global ratio over 1e-7. SimpleFluidProperties has d = 0,
# so it isolates the interface-pressure AD chain (d/dpp_top and d/dsat_n) itself, which is
# exact. (CO2FluidProperties remains correct for actual SOLVES -- the residual is exact and
# a ~1e-3 Jacobian error only mildly slows Newton; see HANDOVER Gotcha D.)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  VEDictator = ve_dictator
[]

[FluidProperties]
  [co2_fp] # CO2-like: light, compressible (bulk_modulus low -> d(rho)/dp non-trivial)
    type = SimpleFluidProperties
    density0 = 700.0
    bulk_modulus = 2.5e7
    viscosity = 5.0e-5
    thermal_expansion = 0.0
  []
  [brine_fp] # brine-like: dense, nearly incompressible
    type = SimpleFluidProperties
    density0 = 1000.0
    bulk_modulus = 2.0e9
    viscosity = 1.0e-3
    thermal_expansion = 0.0
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
    S_wr = 0.1
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '100.0 - 0.1 * x' # sloped top -> nonzero grad(z_T) buoyancy drive
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '80.0 - 0.1 * x' # parallel bottom -> H = 20 m constant
  []
  [pp_top_func]
    type = ParsedFunction
    expression = '12.0e6 - 4.0e4 * x' # 12 MPa -> 8 MPa: supercritical CO2, rho varies
  []
[]

[Variables]
  [pp_top]
  []
  [sat_n]
  []
[]

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_func
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5 # interior: plume present (h > 0) so the interface increment is active
  []
[]

[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      function = z_top_func
    []
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      function = z_bottom_func
    []
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
    type = FunctionDirichletBC
    variable = pp_top
    boundary = left
    function = pp_top_func
  []
  [pp_right]
    type = FunctionDirichletBC
    variable = pp_top
    boundary = right
    function = pp_top_func
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top = z_top
    z_bottom = z_bottom
  []
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
    temperature = 323.15
    eos_reference_depth = interface
    sat_n = sat_n
    z_top = z_top
    z_bottom = z_bottom
    S_wr = 0.1
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

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 1.0
  dt = 1.0
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]
