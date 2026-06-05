# VEFVAdvectiveFlux

!syntax description /FVKernels/VEFVAdvectiveFlux

## Overview

`VEFVAdvectiveFlux` implements the cell-centered FV depth-integrated Darcy advective
flux for one fluid phase:

\begin{equation}
F_c \cdot \hat{n} = -H \, K_{nn} \, \frac{k_{r,c} \, \rho_c}{\mu_c}
  \left( \nabla p_{top} \cdot \hat{n} + \rho_c \, |g| \, \nabla z_T \cdot \hat{n} \right)
\end{equation}

A single kernel serves both the CO2 (`fluid_phase = 0`, `variable = sat_n`) and brine
(`fluid_phase = 1`, `variable = pp_top`) equations.

### Face evaluation via the functor system

All spatially varying quantities are evaluated through MOOSE functors at a
boundary-aware face argument:

- $H = z_{top}(\text{face}) - z_{bottom}(\text{face})$
- $\nabla z_T \cdot \hat{n}$ — `z_top` functor gradient at face
- $\nabla p \cdot \hat{n}$ — `pp_top` functor gradient (carries AD wrt `pp_top`)
- $k_{r,c}$ — `ve_relperm_{n,w}` evaluated at the **upwind** face via `makeFace` (picks
  up Dirichlet BC values on inflow boundaries)
- $\rho_c, \mu_c$ — `ve_density_{n,w}`, `ve_viscosity_{n,w}` evaluated at the
  central-difference (face-averaged) argument, so they are face-correct and carry AD
  wrt `pp_top` on both adjacent cells

### Heterogeneous permeability

The normal permeability $K_{nn}$ uses a distance-weighted harmonic mean of the
element-side and neighbour-side `ve_K_up` tensors on interior faces, giving the
correct TPFA two-point transmissibility across permeability contrasts.  The
element-side value is used on boundary faces.

### Capillary option

When `capillary = true`, the CO2 Darcy potential is augmented by
$\nabla P_c^{up} \cdot \hat{n}$ via the coefficient functors from `VEFVCapPressure`.

Requires `SMP full = true` in `[Preconditioning]` for a coupled FV solve.

## Example Input File Syntax

```text
[FVKernels]
  [co2_flux]
    type = VEFVAdvectiveFlux
    variable   = sat_n
    fluid_phase = 0
    pp_top     = pp_top
    z_top      = z_top
    z_bottom   = z_bottom
  []
  [brine_flux]
    type = VEFVAdvectiveFlux
    variable   = pp_top
    fluid_phase = 1
    pp_top     = pp_top
    z_top      = z_top
    z_bottom   = z_bottom
  []
[]
```

!syntax parameters /FVKernels/VEFVAdvectiveFlux

!syntax inputs /FVKernels/VEFVAdvectiveFlux

!syntax children /FVKernels/VEFVAdvectiveFlux
