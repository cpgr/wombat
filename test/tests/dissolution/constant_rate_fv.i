# Constant-rate dissolution verification -- finite-volume solver.
#
# FV counterpart of constant_rate.i, exercising VEFVDissolutionSink in an FV transient
# solve. CO2 is made immobile (krn_max = 0) so the CO2 flux vanishes identically and the
# CO2 mass equation is one-way decoupled from pressure (sat_n drives pp, not the reverse),
# collapsing to the same pure ODE:
#
#   d/dt(H phi rho_c sat_n) = -q0 * gate(sat_n)
#
# Brine remains mobile (krw_max = 1) and flows in through the open (Dirichlet) pressure
# boundaries to fill the volume freed by dissolution, so the brine equation is consistent.
#
# H=10, phi=0.2, rho_c=700 (H phi rho_c = 1400), q0 = 1.4 kg/m^2/s, sat_n0 = 0.5, s_ref=0.05:
#   sat_n(100 s) = 0.4   (backward-Euler exact: constant RHS)
#   dissolved CO2 mass = q0*t*area = 140 kg
# Conservation: free CO2 lost = H phi rho_c * (0.5-0.4) * area = 140 kg = dissolved.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
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
    krn_max = 0.0 # CO2 immobile: isolates storage + dissolution sink (no advection)
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

[ICs]
  [pp_top_ic]
    type = ConstantIC
    variable = pp_top
    value = 1.0e6
  []
  [sat_n_ic]
    type = ConstantIC
    variable = sat_n
    value = 0.5
  []
[]

[AuxVariables]
  [z_top]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [z_bottom]
    type = MooseVariableFVReal
    initial_condition = -10.0
  []
  [c_diss]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
[]

[AuxKernels]
  [c_diss_aux]
    type = VEDissolvedCO2Aux
    variable = c_diss
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
  [co2_dissolution]
    type = VEFVDissolutionSink
    variable = sat_n
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
  # Open pressure boundaries -> brine reservoir to fill the freed pore volume.
  [pp_left]
    type = FVDirichletBC
    variable = pp_top
    boundary = left
    value = 1.0e6
  []
  [pp_right]
    type = FVDirichletBC
    variable = pp_top
    boundary = right
    value = 1.0e6
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
  [dissolution]
    type = VEDissolution
    dissolution_flux = 1.4
    s_ref = 0.05
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
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 100.0
  dt = 20.0

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  automatic_scaling = true

  nl_rel_tol = 1.0e-9
  nl_abs_tol = 1.0e-11
  nl_max_its = 20
[]

[Postprocessors]
  [sat_n_avg]
    type = ElementAverageValue
    variable = sat_n
    execute_on = 'FINAL'
  []
  [dissolved_co2_mass]
    type = ElementIntegralVariablePostprocessor
    variable = c_diss
    execute_on = 'FINAL'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
