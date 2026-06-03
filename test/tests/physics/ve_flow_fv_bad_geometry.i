# Negative test: FV VE flow physics with z_top wrongly pre-declared as a regular FE
# (LAGRANGE) aux variable. The action declares z_top as MooseVariableFVReal, so the FE
# declaration collides with it and MOOSE rejects the inconsistent variable type -- the
# guard against an FE aux (which is not reinitialised on FV faces, reads zero there, and
# would silently kill the advective flux). See physics/tests (err_fv_geometry_fe_variable).

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
  []
[]

# WRONG on purpose: z_top declared as a continuous FE variable in an FV physics.
[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
  []
[]

[Executioner]
  type = Transient
  dt = 1.0
  num_steps = 1
  solve_type = NEWTON
[]
