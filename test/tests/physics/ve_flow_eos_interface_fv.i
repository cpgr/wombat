# VE flow [Physics] Jacobian -- interface EOS, FV.
#
# Physics-block counterpart of fluidproperties/eos_interface_fv_jacobian.i:
# a fully consistent FV interface-mode model. Both fluid materials evaluate
# density/viscosity at the CO2-brine contact pressure:
#   - VEFVFluidProperties (functor) feeds the FV FLUX (VEFVAdvectiveFlux);
#   - VEFluidProperties (elemental) feeds the FV mass-STORAGE kernel
#     (VEFVMassTimeDerivative, via the stateful ve_density old value).
# PetscJacobianTester confirms d(rho)/d(sat_n) (through h = sat_n*H/(1-S_wr))
# is seeded correctly through both the face evaluation in the flux and the
# heavily-weighted storage term.
#
# SimpleFluidProperties is used for the same reason as the FE interface test:
# interface mode amplifies value-vs-derivative inconsistency by |g|*h. See
# HANDOVER Gotcha D and eos_interface_fe_jacobian.i.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 100
[]

[FluidProperties]
  [co2_fp] # CO2-like: light, compressible
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
    eos_reference_depth = interface
    S_wr = 0.1
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
    value = 0.5 # off the relperm kink, h > 0 so the interface increment is live
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
