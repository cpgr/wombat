# VETrappedSaturationAux

!syntax description /AuxKernels/VETrappedSaturationAux

## Overview

`VETrappedSaturationAux` writes `sat_n_trap`, the depth-averaged residually trapped CO2
saturation, by delegating to the relperm UserObject's `trappedSaturation(sat_n,
sat_n_max)` method.

The result is consumed by `VETrappedCO2MassPostprocessor` and (as a complement)
`VEMobileCO2MassPostprocessor` uses `sat_n - sat_n_trap` for the mobile inventory.

For a non-hysteretic UserObject (`VERelPermSharpUO`, `VERelPermTableUO`),
`trappedSaturation` returns 0, so this aux is a no-op — existing non-hysteretic inputs
are unaffected.

## Example Input File Syntax

```text
[AuxVariables]
  [sat_n_trap]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
[]
[AuxKernels]
  [sat_n_trap_aux]
    type = VETrappedSaturationAux
    variable    = sat_n_trap
    relperm_uo  = relperm_uo
    sat_n       = sat_n
    sat_n_max   = sat_n_max
    execute_on  = TIMESTEP_END
  []
[]
```

!syntax parameters /AuxKernels/VETrappedSaturationAux

!syntax inputs /AuxKernels/VETrappedSaturationAux

!syntax children /AuxKernels/VETrappedSaturationAux
