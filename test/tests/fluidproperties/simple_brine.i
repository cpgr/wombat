# SimpleBrineFluidProperties: fixed-salinity brine through the (p, T) interface.
#
# No flow solve (Problem/solve=false). VEFluidProperties is evaluated at a uniform
# pp_top and temperature, with the wetting phase supplied by SimpleBrineFluidProperties
# (which wraps the three-variable BrineFluidProperties at a constant salt mass fraction)
# and the non-wetting phase by ConstantFluidProperties as an anchor. The brine density
# and viscosity (the ve_density / ve_viscosity wetting-phase components) are extracted to
# auxvariables and reported as postprocessors, regressing the wrapper against the
# underlying BrineFluidProperties EOS.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Problem]
  solve = false
[]

[FluidProperties]
  [co2_fp]
    type = ConstantFluidProperties
    density = 700.0
    viscosity = 5.0e-5
  []
  [brine_fp]
    type = SimpleBrineFluidProperties
    salt_mass_fraction = 0.1
  []
[]

[AuxVariables]
  [pp_top]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 20.0e6
  []
  [co2_rho]
    family = MONOMIAL
    order = CONSTANT
  []
  [co2_mu]
    family = MONOMIAL
    order = CONSTANT
  []
  [brine_rho]
    family = MONOMIAL
    order = CONSTANT
  []
  [brine_mu]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [co2_rho]
    type = ADMaterialStdVectorAux
    variable = co2_rho
    property = ve_density
    index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [co2_mu]
    type = ADMaterialStdVectorAux
    variable = co2_mu
    property = ve_viscosity
    index = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [brine_rho]
    type = ADMaterialStdVectorAux
    variable = brine_rho
    property = ve_density
    index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [brine_mu]
    type = ADMaterialStdVectorAux
    variable = brine_mu
    property = ve_viscosity
    index = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Materials]
  [fluid_props]
    type = VEFluidProperties
    fp_nw = co2_fp
    fp_w = brine_fp
    pp_top = pp_top
    temperature = 323.15
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Postprocessors]
  [co2_density]
    type = ElementAverageValue
    variable = co2_rho
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [co2_viscosity]
    type = ElementAverageValue
    variable = co2_mu
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [brine_density]
    type = ElementAverageValue
    variable = brine_rho
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [brine_viscosity]
    type = ElementAverageValue
    variable = brine_mu
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
