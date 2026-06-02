# Jacobian verification of a fully consistent FV interface-mode model.
#
# Both fluid materials are in eos_reference_depth = interface:
#   - VEFVFluidProperties (functor material) feeds the FV FLUX (VEFVAdvectiveFlux);
#   - VEFluidProperties (elemental material) feeds the FV mass-STORAGE kernel
#     (VEFVMassTimeDerivative, via the stateful ve_density old value).
# Both build H = z_top - z_bottom from the coupled FV geometry variables, so neither
# needs VEGeometry (which FV models do not run). This is the FV counterpart of
# eos_interface_fe_jacobian.i, and -- unlike an earlier version that left storage in
# top_surface mode -- it exercises the interface density in BOTH the flux and the
# heavily-weighted storage term, i.e. a physically consistent model.
#
# The reference pressure is p = pp_top + rho_n(pp_top)*|g|*h with
# h = sat_n*(z_top - z_bottom)/(1-S_wr); density/viscosity carry AD wrt sat_n (through
# h) on top of wrt pp_top. As in the FE test the fluids are SimpleFluidProperties (EXACT
# analytic derivatives): interface mode amplifies an EOS's value-vs-derivative
# inconsistency by |g|*h, so CO2FluidProperties (Newton-inverted) would trip the 1e-7
# storage-weighted threshold -- see eos_interface_fe_jacobian.i and HANDOVER Gotcha D.

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
    expression = '80.0 - 0.1 * x'
  []
  [pp_top_func]
    type = ParsedFunction
    expression = '12.0e6 - 4.0e4 * x' # 12 MPa -> 8 MPa: supercritical CO2, rho varies
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

[ICs]
  [pp_top_ic]
    type = FunctionIC
    variable = pp_top
    function = pp_top_func
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5 # off the sharp-relperm kink, and h > 0 so the interface increment is live
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

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
  []
  [z_bottom]
    type = MooseVariableFVReal
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
    temperature = 323.15
    eos_reference_depth = interface
    sat_n = sat_n
    z_top = z_top
    z_bottom = z_bottom
    S_wr = 0.1
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
