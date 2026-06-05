# Negative test: K_up_xx declared as plain CONSTANT MONOMIAL (not MooseVariableFVReal)
# in an FV VEFlow physics. The VEFlowFV action shall error at add_material with a message
# explaining that FV neighbour reinit requires MooseVariableFVReal.
# See permeability/tests (err_fv_petrophysics_fe_variable).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
  xmin = 0
  xmax = 2
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

[AuxVariables]
  # WRONG on purpose: plain CONSTANT MONOMIAL without type = MooseVariableFVReal.
  # MOOSE does not reinitialise a regular FE aux on FV neighbour elements, so
  # VEPermeability's coupledValue() reads zero on the neighbour side and the
  # harmonic K face average collapses to zero (silently kills the advective flux).
  [K_up_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  # K_up_yy correct so the error isolates to K_up_xx.
  [K_up_yy]
    type = MooseVariableFVReal
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Physics/VEFlow/FiniteVolume]
  [ve]
    z_top = z_top
    z_bottom = z_bottom
    phi_bar = 0.2
    K_up_xx = K_up_xx
    K_up_yy = K_up_yy
    fp_nw = co2_fp
    fp_w = brine_fp
    relperm_uo = relperm_uo
    gravity = '0 0 0'
  []
[]

[Executioner]
  type = Transient
  dt = 1.0
  num_steps = 1
  solve_type = NEWTON
[]
