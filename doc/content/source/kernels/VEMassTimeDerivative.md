# VEMassTimeDerivative

!syntax description /Kernels/VEMassTimeDerivative

## Overview

`VEMassTimeDerivative` implements the depth-integrated mass storage term for one fluid
phase in the FE formulation:

$$
\frac{\partial}{\partial t}\!\left(H \, \bar{\phi} \, \rho_c \, \bar{S}_c\right)
  = \frac{\left(H \, \bar{\phi} \, \rho_c \, \bar{S}_c\right)^{\text{new}}
           - \left(H \, \bar{\phi} \, \rho_c \, \bar{S}_c\right)^{\text{old}}}{\Delta t}
$$

where $H$ = `ve_H`, $\bar{\phi}$ = `ve_phi_bar`, $\rho_c$ = `ve_density[c]`,
$\bar{S}_c$ = `ve_saturation[c]`.

The explicit (new − old) / Δt differencing rather than `_u_dot` ensures the full
nonlinear accumulation $H \phi \rho(p) S(S_n)$ is differenced correctly once a
pressure-dependent EOS is used.

Inheriting from `ADTimeKernelValue` correctly identifies this as a time-derivative
contribution for predictor/corrector schemes and residual-norm reporting.  Jacobian
entries are assembled automatically through MOOSE AD.

One instance is required per fluid phase:

| Phase | `fluid_phase` | `variable` |
|-------|--------------|------------|
| CO2   | 0 | `sat_n` |
| Brine | 1 | `pp_top` |

## Example Input File Syntax

```text
[Kernels]
  [co2_storage]
    type = VEMassTimeDerivative
    variable = sat_n
    fluid_phase = 0
  []
  [brine_storage]
    type = VEMassTimeDerivative
    variable = pp_top
    fluid_phase = 1
  []
[]
```

!syntax parameters /Kernels/VEMassTimeDerivative

!syntax inputs /Kernels/VEMassTimeDerivative

!syntax children /Kernels/VEMassTimeDerivative
