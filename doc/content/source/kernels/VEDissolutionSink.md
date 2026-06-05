# VEDissolutionSink

!syntax description /Kernels/VEDissolutionSink

## Overview

`VEDissolutionSink` removes free-phase CO2 at the convective dissolution rate
`ve_dissolution_rate` (kg/m$^2$/s) computed by `VEDissolution`.  It is assigned to
`variable = sat_n` (the CO2 mass equation) and contributes a positive residual term:

\begin{equation}
R = +\dot{m}_{diss}
\end{equation}

so that together with the storage term $\partial_t(H \phi \rho_c S_n) + \ldots +
\dot{m}_{diss} = 0$ the free CO2 mass decreases (the dissolution is a sink).

The rate is an areal mass flux matching the depth-integrated storage units, so it
enters directly with no $\rho$ scaling.  Jacobian entries (with respect to `sat_n` via
the gate function) are assembled automatically through AD.

This kernel is optional: add it to an input only when dissolution is desired.  Existing
inputs without it are unaffected.  The FV counterpart is `VEFVDissolutionSink`.

Pair with `VEDissolvedCO2Aux` to track the accumulated dissolved inventory and
`VEDissolution` to compute the rate.

## Example Input File Syntax

```text
[Materials]
  [dissolution]
    type = VEDissolution
    dissolution_flux = 1.4
    s_ref = 0.05
  []
[]
[Kernels]
  [co2_dissolution]
    type = VEDissolutionSink
    variable = sat_n
  []
[]
```

!syntax parameters /Kernels/VEDissolutionSink

!syntax inputs /Kernels/VEDissolutionSink

!syntax children /Kernels/VEDissolutionSink
