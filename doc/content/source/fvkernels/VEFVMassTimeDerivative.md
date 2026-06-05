# VEFVMassTimeDerivative

!syntax description /FVKernels/VEFVMassTimeDerivative

## Overview

`VEFVMassTimeDerivative` is the cell-centered FV analogue of `VEMassTimeDerivative`.
It implements the depth-integrated mass storage term for one fluid phase:

\begin{equation}
\frac{\left(H \, \bar{\phi} \, \rho_c \, \bar{S}_c\right)^{\text{new}}
      - \left(H \, \bar{\phi} \, \rho_c \, \bar{S}_c\right)^{\text{old}}}{\Delta t}
\end{equation}

Formation thickness $H$ is built inline from FV variable cell values
$H = z_{top} - z_{bottom}$ rather than from `VEGeometry`'s `ve_H`, because
FE-coupled materials are not reinitialised on FV faces.  This ensures geometry is
sourced consistently with `VEFVAdvectiveFlux`.

One instance per fluid phase:

| Phase | `fluid_phase` | `variable` |
|-------|--------------|------------|
| CO2   | 0 | `sat_n` |
| Brine | 1 | `pp_top` |

Reads `ve_phi_bar` from `VEPorosity`, `ve_density` and `ve_saturation` (with old
values) from `VEFluidProperties` and `VESaturation`.

## Example Input File Syntax

```text
[FVKernels]
  [co2_storage]
    type = VEFVMassTimeDerivative
    variable   = sat_n
    fluid_phase = 0
    z_top      = z_top
    z_bottom   = z_bottom
  []
  [brine_storage]
    type = VEFVMassTimeDerivative
    variable   = pp_top
    fluid_phase = 1
    z_top      = z_top
    z_bottom   = z_bottom
  []
[]
```

!syntax parameters /FVKernels/VEFVMassTimeDerivative

!syntax inputs /FVKernels/VEFVMassTimeDerivative

!syntax children /FVKernels/VEFVMassTimeDerivative
