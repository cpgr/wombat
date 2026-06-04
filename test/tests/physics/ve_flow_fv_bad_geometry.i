# Negative test: FV VE flow physics with define_geometry_variables = false (the user owns
# z_top / z_bottom) but z_top wrongly declared as a regular FE (LAGRANGE) variable. The
# action defers declaration, then checkGeometryVariableType errors because an FE aux is not
# reinitialised on FV faces (reads zero there, would silently kill the advective flux).
# See physics/tests (err_fv_geometry_fe_variable).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
  xmin = 0
  xmax = 1
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

[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = 1.0e-12
    K_up_yy = 1.0e-12
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    define_geometry_variables = false # user provides z_top / z_bottom
  []
[]

[AuxVariables]
  # WRONG on purpose: z_top is a continuous FE variable in an FV physics.
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
  # z_bottom correct (FV) so the error isolates to z_top.
  [z_bottom]
    type = MooseVariableFVReal
  []
[]

[Executioner]
  type = Transient
  dt = 1.0
  num_steps = 1
  solve_type = NEWTON
[]
