# VEAdvectiveFluxP

!syntax description /Kernels/VEAdvectiveFluxP

## Overview

`VEAdvectiveFluxP` implements the depth-integrated Darcy advective flux for the
**pressure variable** (`variable = pp_top`).  It is used on the brine (wetting phase,
`fluid_phase = 1`) mass conservation equation.

After integration by parts the weak-form contribution is:

$$
\int_\Omega \nabla\psi \cdot F_c \, d\Omega
$$

with

$$
F_c = -H \, \mathbf{K}_{up} \, \frac{k_{r,c}^{up}}{\mu_c} \, \rho_c
         \left( \nabla p_{top} + \rho_c \, |g| \, \nabla z_T \right)
$$

Because the kernel variable is `pp_top`, MOOSE's `_grad_u[_qp]` already holds
$\nabla p_{top}$ with full AD derivatives — no extra coupling is needed.

The sloped-surface buoyancy term $\rho_c |g| \nabla z_T$ is the dominant updip
migration mechanism at field scale.  `ve_grad_z_top` is sourced from `VEGeometry`.

See also `VEAdvectiveFluxS` (CO2 equation, on the saturation variable).

## Example Input File Syntax

```text
[Kernels]
  [brine_flux]
    type = VEAdvectiveFluxP
    variable   = pp_top
    fluid_phase = 1
  []
[]
```

!syntax parameters /Kernels/VEAdvectiveFluxP

!syntax inputs /Kernels/VEAdvectiveFluxP

!syntax children /Kernels/VEAdvectiveFluxP
