# VEAdvectiveFluxS

!syntax description /Kernels/VEAdvectiveFluxS

## Overview

`VEAdvectiveFluxS` implements the depth-integrated Darcy advective flux for the
**saturation variable** (`variable = sat_n`).  It is used on the CO2 (non-wetting
phase, `fluid_phase = 0`) mass conservation equation.

Because the kernel variable is `sat_n`, `_grad_u[_qp]` holds $\nabla S_n$ — NOT the
pressure gradient needed for Darcy flow.  The kernel therefore couples explicitly to
`pp_top` via `adCoupledGradient("pp_top")` so off-diagonal Jacobian entries with
respect to pressure are exact.

The flux is:

$$
F_n = -H \, \mathbf{K}_{up} \, \frac{k_{r,n}^{up}}{\mu_n} \, \rho_n
         \left( \nabla p_{top} + \rho_n \, |g| \, \nabla z_T \right)
$$

When `capillary = true` (or `capillary_fringe`), the CO2 Darcy potential is augmented
by $\nabla P_c^{up}$:

$$
\nabla P_c^{up} = \texttt{ve\_dPcup\_dsatn} \cdot \nabla S_n
                  + \texttt{ve\_dPcup\_dH} \cdot \nabla H
$$

sourced from `VEUpscaledCapPressure`.  The `VEPlumeReconstruction` and
`VEUpscaledCapPressure` materials must be present when `capillary != none`.

See also `VEAdvectiveFluxP` (brine equation, on the pressure variable).

## Example Input File Syntax

Advection only:

```text
[Kernels]
  [co2_flux]
    type = VEAdvectiveFluxS
    variable   = sat_n
    fluid_phase = 0
    pp_top      = pp_top
  []
[]
```

With sharp-interface capillary diffusion:

```text
[Kernels]
  [co2_flux]
    type = VEAdvectiveFluxS
    variable   = sat_n
    fluid_phase = 0
    pp_top      = pp_top
    capillary   = true
  []
[]
```

!syntax parameters /Kernels/VEAdvectiveFluxS

!syntax inputs /Kernels/VEAdvectiveFluxS

!syntax children /Kernels/VEAdvectiveFluxS
