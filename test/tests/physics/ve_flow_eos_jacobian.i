# VE flow [Physics] Jacobian -- real EOS (FV, top_surface).
#
# Physics-block counterpart of fluidproperties/fv_eos_jacobian.i: exercises the
# face-correct density path with CO2FluidProperties (non-wetting) and
# SimpleBrineFluidProperties (wetting, fixed salinity). A pressure gradient and a
# sloped top surface make the mass-flux coefficient AND the buoyancy term
# density-dependent. PetscJacobianTester confirms d(rho)/d(pp_top) is seeded
# correctly through the face evaluation in VEFVFluidProperties.
#
# eos_reference_depth = top_surface (default). CO2FluidProperties is safe here:
# the amplification factor |g|*h only matters in interface mode (see
# eos_interface_fv below and HANDOVER Gotcha D).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 100
[]

[FluidProperties]
  [co2_fp]
    type = CO2FluidProperties
  []
  [brine_fp]
    type = SimpleBrineFluidProperties
    salt_mass_fraction = 0.1
  []
[]

[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.1
    krn_max = 1.0
    krw_max = 1.0
  []
[]

[Physics/VEFlow/FiniteVolume]
  [ve]
    define_geometry_variables = false # user declares z_top / z_bottom below
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    temperature = 323.15
    # eos_reference_depth = top_surface (default)
  []
[]

[Functions]
  [z_top_func]
    type = ParsedFunction
    expression = '100.0 - 0.1 * x' # sloped top -> nonzero grad(z_T) buoyancy drive
  []
  [z_bottom_func]
    type = ParsedFunction
    expression = '80.0 - 0.1 * x'
  []
  [pp_top_func]
    type = ParsedFunction
    expression = '12.0e6 - 4.0e4 * x' # 12 MPa -> 8 MPa: supercritical CO2, rho varies
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

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_func
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5 # off the sharp-relperm kink, so CO2 mobility and its derivative are live
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

[FVBCs]
  [pp_left]
    type = FVFunctionDirichletBC
    variable = pp_top
    boundary = left
    function = pp_top_func
  []
  [pp_right]
    type = FVFunctionDirichletBC
    variable = pp_top
    boundary = right
    function = pp_top_func
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
