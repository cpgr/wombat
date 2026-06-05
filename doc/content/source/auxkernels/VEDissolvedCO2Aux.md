# VEDissolvedCO2Aux

!syntax description /AuxKernels/VEDissolvedCO2Aux

## Overview

`VEDissolvedCO2Aux` accumulates the areal dissolved CO2 mass $c_{diss}$ (kg/m$^2$)
removed from the free phase by `VEDissolutionSink`:

\begin{equation}
c_{diss} \leftarrow c_{diss}^{old} + \dot{m}_{diss} \cdot \Delta t
\end{equation}

It executes at `TIMESTEP_END` using the converged end-of-step rate, so the dissolved
inventory gained exactly equals the free CO2 mass removed (mass conservation under
backward-Euler time integration).

Because `c_diss` is frozen within the next solve (lagged/explicit), it can feed the
`VEDissolution` column-capacity cap without adding a nonlinear coupling.

The total dissolved CO2 mass is the domain integral of `c_diss` (use
`ElementIntegralVariablePostprocessor` or `VEDissolvedCO2MassPostprocessor`).

!alert note title=c_diss must be ELEMENTAL
`VEDissolvedCO2Aux` reads `ve_dissolution_rate`, a material property defined at
quadrature points, not nodes.  Declare `c_diss` as `order = CONSTANT, family = MONOMIAL`.

## Example Input File Syntax

```text
[AuxVariables]
  [c_diss]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
[]
[AuxKernels]
  [c_diss_aux]
    type = VEDissolvedCO2Aux
    variable = c_diss
    execute_on = TIMESTEP_END
  []
[]
```

!syntax parameters /AuxKernels/VEDissolvedCO2Aux

!syntax inputs /AuxKernels/VEDissolvedCO2Aux

!syntax children /AuxKernels/VEDissolvedCO2Aux
