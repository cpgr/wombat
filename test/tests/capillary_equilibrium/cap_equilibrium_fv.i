# Capillary no-flow equilibrium on a laterally varying thickness -- FV
#
# FV counterpart of cap_equilibrium_fe.i. See that file for the full rationale.
# VEFVCapPressure declares ve_dPcup_dsatn and ve_dPcup_dH; VEFVAdvectiveFlux with
# capillary=true assembles grad(Pc^up).n = ve_dPcup_dsatn*grad(sat_n).n +
# ve_dPcup_dH*grad(H).n, where grad(H).n comes from the z_top/z_bottom functor
# gradients.
#
#   flat top z_T = 0, H(x) = 10 + 0.2 x  (grad(z_T) = 0, grad(H) = 0.2)
#   equilibrium sat_n_eq(x) = 6 / H(x) in [0.2, 0.6]  ->  sat_n*H = const
#   grad(Pc^up) = delta_rho*g*grad(sat_n*H) = 0  ->  no-flow fixed point.
#
# satn_drift (ElementL2Error vs sat_n_eq) stays ~O(1e-4) WITH the grad(H) term
# and grows to O(1e-2)+ without it.

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

[Functions]
  [z_bottom_fn]
    type = ParsedFunction
    expression = '-(10.0 + 0.2 * x)'
  []
  [sat_n_eq]
    type = ParsedFunction
    expression = '6.0 / (10.0 + 0.2 * x)'
  []
[]

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e5
  []
  [sat_n_ic]
    type = FunctionIC
    variable = sat_n
    function = sat_n_eq
  []
  [z_bottom_ic]
    type = FunctionIC
    variable = z_bottom
    function = z_bottom_fn
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
    initial_condition = 0.0
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

[FVBCs]
  # Single pressure pin removes the closed-domain pressure null space.
  [pp_pin]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e5
  []
[]

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

[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n = sat_n
  []
  [cap_pressure]
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

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 2.0e4
  dt = 2.0e3

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
  [satn_drift]
    type = ElementL2Error
    variable = sat_n
    function = sat_n_eq
    execute_on = 'FINAL'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
