# VEGeometry

!syntax description /Materials/VEGeometry

## Overview

`VEGeometry` computes the three geometric material properties consumed by all other VE
materials and kernels:

- `ve_H` — formation thickness $H = z_T - z_B$ (m)
- `ve_grad_z_top` — 2-D horizontal gradient of the top-surface elevation $\nabla z_T$ (m/m)
- `ve_grad_H` — 2-D horizontal gradient of the formation thickness $\nabla H$ (m/m)

`ve_grad_z_top` enters every `VEAdvectiveFluxS`/`VEAdvectiveFluxP` kernel as the
sloped-surface buoyancy drive $\rho_c \, g \, \nabla z_T$ — the dominant updip CO2
migration mechanism at field scale.  `ve_grad_H` provides the $\nabla H$ term for the
chain-rule gradient of the upscaled capillary pressure when `capillary = true`.

All three properties are plain `Real` (not `ADMaterialProperty`) because formation
geometry is fixed in time and carries no Jacobian contribution with respect to the
primary flow variables.

### Geometry variables

Both `z_top` and `z_bottom` must be `LAGRANGE FIRST ORDER` `AuxVariables`.
`VEGeometry` calls `coupledGradient` to compute $\nabla z_T$ and $\nabla z_B$; a
`CONSTANT MONOMIAL` declaration would produce a zero gradient and silently kill the
entire buoyancy drive.

The `ex2ve` preprocessing tool emits nodal `z_top`/`z_bottom` fields directly from
the Petrel-to-2D upscaling workflow.  Formation thickness is computed internally as
$H = z_T - z_B$; a precomputed thickness field is **not** accepted.

!alert note title=FV models do not use VEGeometry
The FV kernels (`VEFVAdvectiveFlux`, `VEFVMassTimeDerivative`) build $H$ and
$\nabla z_T$ inline from the `z_top`/`z_bottom` FV variable values.  `VEGeometry` is
only required in FE (`VEAdvectiveFluxS`, `VEAdvectiveFluxP`, `VEMassTimeDerivative`)
inputs.  The `VEFlowFE` Physics action adds `VEGeometry` automatically; the `VEFlowFV`
action does not.

## Example Input File Syntax

```text
[AuxVariables]
  [z_top]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0   # or read from Exodus
  []
  [z_bottom]
    order = FIRST
    family = LAGRANGE
    initial_condition = -50.0
  []
[]

[Materials]
  [geometry]
    type = VEGeometry
    z_top    = z_top
    z_bottom = z_bottom
  []
[]
```

!syntax parameters /Materials/VEGeometry

!syntax inputs /Materials/VEGeometry

!syntax children /Materials/VEGeometry
