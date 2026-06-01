# Analytic verification of the Land residual-trapping algebra and the sat_n_max /
# sat_n_trap state machinery, over a drainage -> imbibition cycle.
#
# No flow solve (Problem/solve=false): sat_n is driven by a prescribed triangle in
# time and the hysteresis bookkeeping is checked against closed-form Land values.
# sat_n rises 0 -> 0.5 (drainage) then falls 0.5 -> 0.2 (imbibition).
#
# With S_wr = 0 (so s = sat_n) and S_gr_max = 0.25, the Land coefficient is
#   C = 1/S_gr_max - 1 = 3,
# and the turning point at sat_n_max = 0.5 gives the residual
#   s_gr = s_max / (1 + C s_max) = 0.5 / (1 + 3*0.5) = 0.2.
# On the imbibition branch the trapped saturation is
#   sat_n_trap = s_gr (s_max - sat_n)/(s_max - s_gr) = 0.2 (0.5 - sat_n)/0.3,
# which is 0 at the turning point and grows to s_gr = 0.2 (all CO2 immobile) when
# sat_n falls to s_gr. Checkpoints (exact):
#   sat_n = 0.5 -> trap 0    (turning point, drainage->imbibition)
#   sat_n = 0.35 -> trap 0.1
#   sat_n = 0.2  -> trap 0.2 (fully trapped; mobile = sat_n - trap = 0)

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  solve = false
[]

[UserObjects]
  [relperm_uo]
    type = VERelPermHysteresisUO
    S_wr = 0.0
    krn_max = 1.0
    krw_max = 1.0
    S_gr_max = 0.25
  []
[]

[Functions]
  [sat_n_drive]
    type = ParsedFunction
    # rise to 0.5 (drainage, t in [0,5]) then fall by 0.15/step (imbibition)
    expression = 'if(t <= 5, 0.1 * t, 0.5 - 0.15 * (t - 5))'
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
  [sat_n_trap]
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
  [sat_n_trap_ic]
    type = ConstantIC
    variable = sat_n_trap
    value = 0.0
  []
[]

[AuxKernels]
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
  [sat_n_trap_aux]
    type = VETrappedSaturationAux
    variable = sat_n_trap
    relperm_uo = relperm_uo
    sat_n = sat_n
    sat_n_max = sat_n_max
    execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 7.0
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
  [sat_n_trap]
    type = ElementAverageValue
    variable = sat_n_trap
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
