# VERelPerm

!syntax description /Materials/VERelPerm

## Overview

`VERelPerm` is the FE adapter that evaluates the upscaled relative permeabilities
from a `VERelativePermeability` UserObject and publishes the result as the
`ADMaterialProperty<std::vector<Real>>` named `ve_relperm`:

- `ve_relperm[0]` — CO2 (non-wetting) relative permeability
- `ve_relperm[1]` — brine (wetting) relative permeability

It reads `ve_saturation` from `VESaturation` and passes the current `sat_n` value
(with AD derivative seeding) to the UserObject's `relativePermeabilityAD` method.
Swapping the UserObject (e.g. `VERelPermSharpUO` → `VERelPermTableUO` →
`VERelPermHysteresisUO`) changes the kr model without touching any kernel.

For hysteretic models, optionally couple the `sat_n_max` AuxVariable written by
`VESaturationMaxAux`.  When provided, it is passed as a frozen (plain `Real`) turning
point to the three-argument scanning-curve overload.  Existing non-hysteretic inputs
are unaffected.

The FV counterpart is `VEFVRelPerm`.

## Example Input File Syntax

Non-hysteretic (verification):

```text
[UserObjects]
  [relperm_uo]
    type = VERelPermSharpUO
    S_wr = 0.2
    krn_max = 1.0
    krw_max = 1.0
  []
[]
[Materials]
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
  []
[]
```

Hysteretic:

```text
[Materials]
  [relperm]
    type = VERelPerm
    relperm_uo = relperm_uo
    sat_n_max  = sat_n_max   # AuxVariable from VESaturationMaxAux
  []
[]
```

!syntax parameters /Materials/VERelPerm

!syntax inputs /Materials/VERelPerm

!syntax children /Materials/VERelPerm
