# VEFVRelPerm

!syntax description /Materials/VEFVRelPerm

## Overview

`VEFVRelPerm` is the FV adapter: a `FunctorMaterial` that publishes the upscaled
relative permeabilities as the functors

- `ve_relperm_n` — CO2 relative permeability
- `ve_relperm_w` — brine relative permeability

delegating to a `VERelativePermeability` UserObject.  The functors evaluate `sat_n` at
whatever face argument `VEFVAdvectiveFlux` supplies (typically the upwind cell for the
advected quantity), so a Dirichlet `sat_n` inlet BC drives inflow correctly.

For hysteretic models, supply the optional `sat_n_max` functor (written by
`VESaturationMaxAux`): it is evaluated at the same face and passed as a frozen `Real`
turning point to the scanning-curve overload.  Existing inputs without `sat_n_max` are
unaffected.

The FE counterpart is `VERelPerm`.

## Example Input File Syntax

```text
[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n      = sat_n
  []
[]
```

With hysteresis:

```text
[FunctorMaterials]
  [relperm]
    type = VEFVRelPerm
    relperm_uo = relperm_uo
    sat_n      = sat_n
    sat_n_max  = sat_n_max
  []
[]
```

!syntax parameters /Materials/VEFVRelPerm

!syntax inputs /Materials/VEFVRelPerm

!syntax children /Materials/VEFVRelPerm
