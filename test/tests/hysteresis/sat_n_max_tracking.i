# Phase 1 diagnostic: VESaturationMaxAux tracks the running maximum of sat_n.
#
# No flow solve (Problem/solve=false). sat_n is driven by a prescribed triangle
# profile in time -- it rises 0 -> 0.6 (t in [0,5], drainage) then falls 0.6 ->
# 0.2 (t in [5,10], imbibition). VESaturationMaxAux must hold the running max:
# sat_n_max climbs with sat_n on the way up, then stays pinned at the 0.6 peak as
# sat_n recedes. This is the turning point the hysteretic relperm will use; here
# it is purely diagnostic (nothing consumes it yet).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  solve = false
[]

[Functions]
  [sat_n_drive]
    type = ParsedFunction
    # rise (drainage) then fall (imbibition); peak 0.6 at t = 5
    expression = 'if(t <= 5, 0.12 * t, 0.6 - 0.08 * (t - 5))'
  []
[]

[AuxVariables]
  [sat_n]
    family = MONOMIAL
    order = CONSTANT
  []
  [sat_n_max]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [sat_n_max_ic]
    type = ConstantIC
    variable = sat_n_max
    value = 0.0
  []
[]

[AuxKernels]
  # Update sat_n at the start of each step so VESaturationMaxAux (TIMESTEP_END)
  # sees the current value.
  [sat_n_drive_aux]
    type = FunctionAux
    variable = sat_n
    function = sat_n_drive
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [sat_n_max_aux]
    type = VESaturationMaxAux
    variable = sat_n_max
    sat_n = sat_n
    # execute_on defaults to TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 10.0
  dt = 1.0
[]

[Postprocessors]
  [sat_n]
    type = ElementAverageValue
    variable = sat_n
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sat_n_max]
    type = ElementAverageValue
    variable = sat_n_max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
