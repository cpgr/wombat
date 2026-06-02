# VEFVAdvectiveOutflowBC verification -- a phase leaves the domain instead of piling
# up against a no-flow wall.
#
# 1D tilted column with CO2 placed ONLY in the cell against the open updip (+x)
# boundary (sat_n = 0.5 for x > 98, else 0). That boundary carries a Dirichlet brine
# pressure plus VEFVAdvectiveOutflowBC on sat_n; buoyancy (rho_n * g * grad(z_T)) drives
# the CO2 out through it. Because there is no CO2 upstream to refill the cell, it drains
# directly:
#   - co2_integral (domain CO2 volume) DECREASES monotonically from its initial
#     0.5 * 2 = 1.0, and
#   - sat_at_outflow (that same boundary cell) falls steadily from 0.5,
# as CO2 leaves through the outflow BC. A no-flow sat_n BC would hold both FIXED -- that
# contrast is the verification signature. The drain is gradual, not a fast empty-out:
# with the downdip end closed it is paced by counter-current brine inflow at the single
# open boundary and by the interior kr_n carried in the outflow flux. Over this 1e6 s
# window sat_at_outflow declines ~0.5 -> ~0.40. The exact rate is set by the coupled
# nonlinear flux and is mildly timestep-sensitive, so this is a REGRESSION check on
# monotone outflow, not an analytical decay rate -- regenerate the gold if dt changes.
#
# Why only the boundary cell (not a wider blob): the downdip end is CLOSED, so CO2 can
# only leave by slow counter-current exchange at the single open boundary (brine must
# re-enter as CO2 exits). A wider blob keeps the boundary cell topped up from upstream
# far longer than a short test runs, so sat_at_outflow would sit at 0.5. Isolating the
# CO2 at the boundary shows the drain directly. The signature checked is that, with the
# outflow BC, both quantities fall; the default natural (no-flow) BC would conserve them
# (CO2 piling up at x = 100).
#
# Parameters mirror nordbotten_celia/nc_sloped_fv (K=1e-12, phi=0.2, H=20, theta=2deg).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 50
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  VEDictator = ve_dictator
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
  # CO2 only in the cell against the outflow boundary (centroid x = 99, dx = 2).
  [sat_n_init]
    type = ParsedFunction
    expression = 'if(x > 98, 0.5, 0.0)'
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
    function = pp_top_init
  []
  [sat_n_ic]
    type = FunctionIC
    variable = sat_n
    function = sat_n_init
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
  # Open updip (right) boundary: brine pressure Dirichlet + CO2 advective outflow.
  [pp_right]
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
    z_top = z_top
    z_bottom = z_bottom
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
    type = VEFluidPropertiesConst
    rho_co2 = 700.0
    rho_brine = 1020.0
    mu_co2 = 5.0e-5
    mu_brine = 8.0e-4
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
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  # 20 steps over 1e6 s: a clear monotone drain (sat_at_outflow ~0.5 -> ~0.40) without a
  # long flat tail. The drain is gradual (counter-current limited); extend end_time for a
  # larger drop. The rate is mildly timestep-sensitive, so regenerate the gold if dt
  # changes.
  end_time = 1.0e6
  dt = 5.0e4

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8
  nl_max_its = 25
  l_tol = 1.0e-6
[]

[Postprocessors]
  # Domain CO2 volume. Starts at 0.5 * 2 = 1.0; must DECREASE as CO2 exits the updip
  # outflow boundary. A no-flow sat_n BC would hold it at 1.0.
  [co2_integral]
    type = ElementIntegralVariablePostprocessor
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # CO2 saturation in the boundary cell itself. Starts at 0.5 and declines steadily as
  # it drains out (no upstream CO2 to refill it).
  [sat_at_outflow]
    type = PointValue
    variable = sat_n
    point = '99 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
