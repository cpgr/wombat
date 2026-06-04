# VE flow [Physics] Jacobian -- interface EOS, FE.
#
# Physics-block counterpart of fluidproperties/eos_interface_fe_jacobian.i:
# eos_reference_depth = interface makes VEFluidProperties evaluate the EOS at
# the CO2-brine contact z_T + h, with the hydrostatic pressure
# p = pp_top + rho_n(pp_top)*|g|*h and h = sat_n*H/(1-S_wr). This adds a
# d(rho)/d(sat_n) coupling (through h) that top_surface mode does not have.
# PetscJacobianTester confirms both d/d(pp_top) and d/d(sat_n) are seeded
# correctly into the storage and flux kernels.
#
# SimpleFluidProperties (exact analytic d(rho)/dp) is used instead of
# CO2FluidProperties. Interface mode amplifies an EOS's value-vs-derivative
# inconsistency by |g|*h (~100 here); CO2's Newton-inverted Span-Wagner has a
# tiny inconsistency that passes top_surface tests but trips the 1e-7 ratio
# threshold in interface mode. SimpleFluidProperties has d = 0 and isolates
# the interface-pressure AD chain. See HANDOVER Gotcha D.

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

# The Physics action creates pp_top, sat_n (FE LAGRANGE) and z_top, z_bottom
# (FE LAGRANGE, define_geometry_variables = true default), plus the full FE
# material chain with eos_reference_depth = interface.
[Physics/VEFlow/FiniteElement]
  [ve]
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
    expression = '80.0 - 0.1 * x' # parallel bottom -> H = 20 m constant
  []
  [pp_top_func]
    type = ParsedFunction
    expression = '12.0e6 - 4.0e4 * x' # 12 MPa -> 8 MPa: supercritical CO2, rho varies
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
  # z_top / z_bottom declared by the Physics action (FE LAGRANGE).
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
