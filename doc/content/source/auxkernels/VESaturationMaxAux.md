# VESaturationMaxAux

!syntax description /AuxKernels/VESaturationMaxAux

## Overview

`VESaturationMaxAux` tracks the maximum CO2 saturation each point has ever experienced —
the turning point of the imbibition scanning curve used by `VERelPermHysteresisUO` for
residual trapping:

$$
S_{n,\max} \leftarrow \max(S_n, S_{n,\max}^{\text{old}})
$$

It is advanced once per timestep (`execute_on = TIMESTEP_END`, after the nonlinear
solve converges), keeping `sat_n_max` frozen within each solve.  This makes the
hysteresis state smooth and explicit — a common, numerically stable design for
scanning-curve hysteresis.

The auxiliary variable must have the same family as `sat_n` (LAGRANGE for FE,
MONOMIAL / `MooseVariableFVReal` for FV) and its initial condition must equal the
`sat_n` IC.

`sat_n_max` is consumed by:

- `VERelPerm` / `VEFVRelPerm` (when `sat_n_max` is coupled) — history-aware kr
- `VETrappedSaturationAux` — computes `sat_n_trap` for the mass postprocessors

## Example Input File Syntax

```text
[AuxVariables]
  [sat_n_max]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0   # should match sat_n IC
  []
[]
[AuxKernels]
  [sat_n_max_aux]
    type = VESaturationMaxAux
    variable    = sat_n_max
    sat_n       = sat_n
    execute_on  = TIMESTEP_END
  []
[]
```

!syntax parameters /AuxKernels/VESaturationMaxAux

!syntax inputs /AuxKernels/VESaturationMaxAux

!syntax children /AuxKernels/VESaturationMaxAux
