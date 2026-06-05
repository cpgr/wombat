# VEFVDissolutionSink

!syntax description /FVKernels/VEFVDissolutionSink

## Overview

`VEFVDissolutionSink` is the FV analogue of `VEDissolutionSink`.  It removes free-phase
CO2 at the convective dissolution rate `ve_dissolution_rate` (kg/m$^2$/s) from
`VEDissolution`, assigned to `variable = sat_n` (the CO2 mass equation).

The residual contribution is $+\dot{m}_{diss}$ (same sign convention as
`VEFVMassTimeDerivative`).  The rate is an areal mass flux matching the
depth-integrated storage units, so it enters directly with no $\rho$ scaling.
Jacobian entries (wrt `sat_n` via the gate) are assembled via AD.

Pair with `VEDissolvedCO2Aux` to track the dissolved inventory.

## Example Input File Syntax

```text
[FVKernels]
  [co2_dissolution]
    type = VEFVDissolutionSink
    variable = sat_n
  []
[]
```

!syntax parameters /FVKernels/VEFVDissolutionSink

!syntax inputs /FVKernels/VEFVDissolutionSink

!syntax children /FVKernels/VEFVDissolutionSink
